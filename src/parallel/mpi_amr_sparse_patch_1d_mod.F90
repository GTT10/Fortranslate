module mpi_amr_sparse_patch_1d_mod
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use transport_database_mod, only: gas_transport_species
  use simulation_config_reactive_1d_mod, only: reactive_1d_config
  use reactive_1d_mod, only: advance_reactive_chemistry
  use amr_hierarchy_1d_mod, only: &
    amr_two_level_hierarchy_1d, amr_level_field_1d, &
    accumulate_coarse_flux_1d, &
    accumulate_fine_flux_1d
  use amr_multipatch_1d_mod, only: &
    prolong_patch_set_1d, average_down_patch_set_1d, &
    synchronize_patch_set_1d
  use amr_regrid_1d_mod, only: amr_regrid_plan_collection_1d
  use amr_reactive_1d_mod, only: &
    recover_level_temperatures_1d, fill_physical_ghosts_1d, &
    advance_amr_level_1d, advance_transport_level_1d
  use amr_patch_tree_1d_mod, only: &
    amr_child_patch_plan_1d, amr_patch_level_plan_1d, &
    amr_patch_tree_hierarchy_1d, &
    amr_patch_tree_relation_flux_registers_1d, &
    initialize_patch_tree_1d, initialize_patch_tree_flux_registers_1d, &
    patch_tree_child_geometry_1d
  use amr_patch_tree_reactive_1d_mod, only: &
    amr_patch_tree_reactive_patch_1d, &
    amr_patch_tree_reactive_solution_1d, fill_one_child_ghosts, &
    plan_tagged_reactive_parent_1d
  use mpi_amr_patch_1d_mod, only: &
    mpi_amr_patch_distribution_1d, &
    initialize_mpi_amr_patch_distribution_1d, &
    mpi_amr_distribution_matches_hierarchy_1d, &
    synchronize_owned_patch_tree_reactive_1d
  implicit none
  private

  integer, parameter :: sparse_patch_migration_tag = 2601
  integer, parameter :: sparse_adjacent_halo_tag = 2602
  integer, parameter :: sparse_child_parent_tag = 2603
  integer, parameter :: sparse_parent_state_tag = 2604
  integer, parameter :: sparse_interval_state_tag = 2605
  integer, parameter :: sparse_boundary_flux_tag = 2606
  integer, parameter :: sparse_flux_correction_tag = 2607
  integer, parameter :: sparse_regrid_prolongation_tag = 2608
  integer, parameter :: sparse_regrid_overlap_tag = 2609

  type, public :: mpi_amr_sparse_communication_counts_1d
    integer :: interval_state_transfers = 0
    integer :: boundary_flux_transfers = 0
    integer :: shared_flux_correction_transfers = 0
  end type mpi_amr_sparse_communication_counts_1d

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
  public :: regrid_sparse_patch_tree_reactive_1d
  public :: regrid_tagged_sparse_patch_tree_reactive_1d
  public :: advance_sparse_patch_tree_chemistry_1d
  public :: advance_sparse_patch_tree_hydro_1d
  public :: advance_sparse_patch_tree_transport_1d
  public :: advance_sparse_patch_tree_reactive_1d

contains

  subroutine advance_sparse_patch_tree_chemistry_1d( &
      species, reactions, config, interval, distribution, solution, ok, &
      local_patch_advances, local_halo_transfers, local_parent_transfers, &
      local_parent_state_transfers)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(reactive_1d_config), intent(in) :: config
    real(dp), intent(in) :: interval
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_patch_advances
    integer, intent(out), optional :: local_halo_transfers
    integer, intent(out), optional :: local_parent_transfers
    integer, intent(out), optional :: local_parent_state_transfers

    type(mpi_amr_sparse_reactive_solution_1d) :: backup
    character(len=32) :: boundary
    logical :: local_ok, accepted, mpi_ok
    integer :: level, patch, nx, advances, halo_transfers, parent_transfers
    integer :: parent_state_transfers

    ok = .false.
    advances = 0
    halo_transfers = 0
    parent_transfers = 0
    parent_state_transfers = 0
    if (present(local_patch_advances)) local_patch_advances = 0
    if (present(local_halo_transfers)) local_halo_transfers = 0
    if (present(local_parent_transfers)) local_parent_transfers = 0
    if (present(local_parent_state_transfers)) &
      local_parent_state_transfers = 0
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
      species, distribution, solution, local_ok, parent_transfers)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) then
      solution = backup
      advances = 0
      return
    end if
    call refresh_sparse_reactive_ghosts_1d( &
      species, config, distribution, solution, local_ok, halo_transfers, &
      parent_state_transfers)
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
    if (present(local_halo_transfers)) &
      local_halo_transfers = halo_transfers
    if (present(local_parent_transfers)) &
      local_parent_transfers = parent_transfers
    if (present(local_parent_state_transfers)) &
      local_parent_state_transfers = parent_state_transfers
  end subroutine advance_sparse_patch_tree_chemistry_1d

  subroutine advance_sparse_patch_tree_hydro_1d( &
      species, config, dt, distribution, solution, ok, &
      local_patch_advances, local_communication)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    real(dp), intent(in) :: dt
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_patch_advances
    type(mpi_amr_sparse_communication_counts_1d), intent(out), optional :: &
      local_communication

    type(mpi_amr_sparse_reactive_solution_1d) :: backup
    type(mpi_amr_sparse_communication_counts_1d) :: communication
    type(amr_patch_tree_relation_flux_registers_1d), allocatable :: registers(:)
    real(dp), allocatable :: left_integral(:), right_integral(:)
    logical :: local_ok, accepted, mpi_ok
    integer :: advances

    ok = .false.
    advances = 0
    if (present(local_patch_advances)) local_patch_advances = 0
    if (present(local_communication)) &
      local_communication = mpi_amr_sparse_communication_counts_1d()
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
      left_integral, right_integral, distribution%owner_of(0, 1), advances, &
      communication, local_ok)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) then
      solution = backup
      advances = 0
      return
    end if
    call synchronize_sparse_counter_delta_1d( &
      distribution, backup%level_advances, solution%level_advances, local_ok)
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
    if (present(local_communication)) local_communication = communication
  end subroutine advance_sparse_patch_tree_hydro_1d

  recursive subroutine advance_sparse_patch_hydro_recursive_1d( &
      species, config, distribution, solution, registers, level, &
      parent_patch, interval, left_integral, right_integral, &
      result_owner, local_advances, communication, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(inout) :: solution
    type(amr_patch_tree_relation_flux_registers_1d), &
      intent(inout) :: registers(:)
    integer, intent(in) :: level, parent_patch
    real(dp), intent(in) :: interval
    real(dp), intent(out) :: left_integral(:), right_integral(:)
    integer, intent(in) :: result_owner
    integer, intent(inout) :: local_advances
    type(mpi_amr_sparse_communication_counts_1d), intent(inout) :: &
      communication
    logical, intent(out) :: ok

    real(dp), allocatable :: state_start(:, :), state_end(:, :), flux(:, :)
    real(dp), allocatable :: child_left(:, :), child_right(:, :)
    real(dp) :: child_interval, alpha, dx
    logical :: local_ok, patch_ok, accepted, mpi_ok, physical_boundary
    integer :: nx, owner, ratio, substep, child, child_count, global_child
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
    patch_ok = .true.
    if (distribution%rank == owner) then
      allocate(state_start(solution%nvar, 0:nx + 1))
      allocate(state_end(solution%nvar, 0:nx + 1))
      allocate(flux(solution%nvar, 0:nx))
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
        left_integral = interval * flux(:, 0)
        right_integral = interval * flux(:, nx)
      end if
    end if
    call all_ranks_accept_sparse_1d( &
      distribution, patch_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) return
    call transfer_sparse_boundary_flux_integrals_1d( &
      distribution, owner, result_owner, left_integral, right_integral, &
      communication, local_ok)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) return
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
    local_ok = .true.
    if (distribution%rank == owner) then
      do child = 1, child_count
        call accumulate_coarse_flux_1d( &
          registers(level)%parents(parent_patch)%children(child), &
          flux(:, solution%hierarchy%relations(level)% &
            child_sets(parent_patch)%patches(child)%fine_coarse_lower - 1), &
          flux(:, solution%hierarchy%relations(level)% &
            child_sets(parent_patch)%patches(child)%fine_coarse_upper), &
          interval, local_ok)
        if (.not. local_ok) exit
      end do
    end if
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) return
    call distribute_sparse_interval_states_to_child_owners_1d( &
      distribution, solution, level, parent_patch, owner, state_start, &
      state_end, communication, local_ok)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) return
    allocate(child_left(solution%nvar, child_count))
    allocate(child_right(solution%nvar, child_count))
    child_left = 0.0_dp
    child_right = 0.0_dp

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
          child_right(:, child), owner, local_advances, communication, &
          local_ok)
        if (.not. local_ok) return
      end do
      call reconcile_sparse_adjacent_child_fluxes_1d( &
        species, distribution, solution, level, parent_patch, child_left, &
        child_right, communication, local_ok)
      call all_ranks_accept_sparse_1d( &
        distribution, local_ok, accepted, mpi_ok)
      if (.not. mpi_ok .or. .not. accepted) return
      local_ok = .true.
      if (distribution%rank == owner) then
        do child = 1, child_count
          call accumulate_fine_flux_1d( &
            registers(level)%parents(parent_patch)%children(child), &
            child_left(:, child) / child_interval, &
            child_right(:, child) / child_interval, child_interval, local_ok)
          if (.not. local_ok) exit
        end do
      end if
      call all_ranks_accept_sparse_1d( &
        distribution, local_ok, accepted, mpi_ok)
      if (.not. mpi_ok .or. .not. accepted) return
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
      local_patch_advances, local_communication)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(reactive_1d_config), intent(in) :: config
    real(dp), intent(in) :: interval
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_patch_advances
    type(mpi_amr_sparse_communication_counts_1d), intent(out), optional :: &
      local_communication

    type(mpi_amr_sparse_reactive_solution_1d) :: backup
    type(mpi_amr_sparse_communication_counts_1d) :: communication
    type(amr_patch_tree_relation_flux_registers_1d), allocatable :: registers(:)
    real(dp), allocatable :: left_integral(:), right_integral(:)
    logical :: local_ok, accepted, mpi_ok
    integer :: advances

    ok = .false.
    advances = 0
    if (present(local_patch_advances)) local_patch_advances = 0
    if (present(local_communication)) &
      local_communication = mpi_amr_sparse_communication_counts_1d()
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
      interval, left_integral, right_integral, distribution%owner_of(0, 1), &
      advances, communication, local_ok)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) then
      solution = backup
      advances = 0
      return
    end if
    call synchronize_sparse_counter_delta_1d( &
      distribution, backup%transport_level_advances, &
      solution%transport_level_advances, local_ok)
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
    if (present(local_communication)) local_communication = communication
  end subroutine advance_sparse_patch_tree_transport_1d

  subroutine advance_sparse_patch_tree_reactive_1d( &
      species, reactions, config, dt, distribution, solution, ok, transport, &
      local_chemistry_advances, local_hydro_advances, &
      local_transport_advances)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(reactive_1d_config), intent(in) :: config
    real(dp), intent(in) :: dt
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok
    type(gas_transport_species), intent(in), optional :: transport(:)
    integer, intent(out), optional :: local_chemistry_advances
    integer, intent(out), optional :: local_hydro_advances
    integer, intent(out), optional :: local_transport_advances

    type(mpi_amr_sparse_reactive_solution_1d) :: backup
    logical :: local_ok, accepted, mpi_ok
    integer :: chemistry_advances, hydro_advances, transport_advances
    integer :: stage_advances

    ok = .false.
    chemistry_advances = 0
    hydro_advances = 0
    transport_advances = 0
    if (present(local_chemistry_advances)) local_chemistry_advances = 0
    if (present(local_hydro_advances)) local_hydro_advances = 0
    if (present(local_transport_advances)) local_transport_advances = 0
    local_ok = dt > 0.0_dp .and. size(species) >= 1 .and. &
      solution%is_valid(distribution)
    if (config%transport_enabled) then
      local_ok = local_ok .and. present(transport)
      if (present(transport)) &
        local_ok = local_ok .and. size(transport) == size(species)
    end if
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) return
    backup = solution

    if (config%chemistry_enabled) then
      call advance_sparse_patch_tree_chemistry_1d( &
        species, reactions, config, 0.5_dp * dt, distribution, solution, &
        local_ok, stage_advances)
      call all_ranks_accept_sparse_1d( &
        distribution, local_ok, accepted, mpi_ok)
      if (.not. mpi_ok .or. .not. accepted) go to 900
      chemistry_advances = chemistry_advances + stage_advances
    end if
    if (config%transport_enabled) then
      call advance_sparse_patch_tree_transport_1d( &
        species, transport, config, 0.5_dp * dt, distribution, solution, &
        local_ok, stage_advances)
      call all_ranks_accept_sparse_1d( &
        distribution, local_ok, accepted, mpi_ok)
      if (.not. mpi_ok .or. .not. accepted) go to 900
      transport_advances = transport_advances + stage_advances
    end if
    call advance_sparse_patch_tree_hydro_1d( &
      species, config, dt, distribution, solution, local_ok, stage_advances)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) go to 900
    hydro_advances = hydro_advances + stage_advances
    if (config%transport_enabled) then
      call advance_sparse_patch_tree_transport_1d( &
        species, transport, config, 0.5_dp * dt, distribution, solution, &
        local_ok, stage_advances)
      call all_ranks_accept_sparse_1d( &
        distribution, local_ok, accepted, mpi_ok)
      if (.not. mpi_ok .or. .not. accepted) go to 900
      transport_advances = transport_advances + stage_advances
    end if
    if (config%chemistry_enabled) then
      call advance_sparse_patch_tree_chemistry_1d( &
        species, reactions, config, 0.5_dp * dt, distribution, solution, &
        local_ok, stage_advances)
      call all_ranks_accept_sparse_1d( &
        distribution, local_ok, accepted, mpi_ok)
      if (.not. mpi_ok .or. .not. accepted) go to 900
      chemistry_advances = chemistry_advances + stage_advances
    end if
    call refresh_sparse_reactive_ghosts_1d( &
      species, config, distribution, solution, local_ok)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) go to 900
    local_ok = solution%is_valid(distribution)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) go to 900
    if (present(local_chemistry_advances)) &
      local_chemistry_advances = chemistry_advances
    if (present(local_hydro_advances)) &
      local_hydro_advances = hydro_advances
    if (present(local_transport_advances)) &
      local_transport_advances = transport_advances
    ok = .true.
    return

900 continue
    solution = backup
    if (present(local_chemistry_advances)) local_chemistry_advances = 0
    if (present(local_hydro_advances)) local_hydro_advances = 0
    if (present(local_transport_advances)) local_transport_advances = 0
    ok = .false.
  end subroutine advance_sparse_patch_tree_reactive_1d

  recursive subroutine advance_sparse_patch_transport_recursive_1d( &
      species, transport, config, distribution, solution, registers, level, &
      parent_patch, interval, left_integral, right_integral, &
      result_owner, local_advances, communication, ok)
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
    integer, intent(in) :: result_owner
    integer, intent(inout) :: local_advances
    type(mpi_amr_sparse_communication_counts_1d), intent(inout) :: &
      communication
    logical, intent(out) :: ok

    real(dp), allocatable :: state_start(:, :), state_end(:, :), flux(:, :)
    real(dp), allocatable :: child_left(:, :), child_right(:, :)
    real(dp) :: child_interval, alpha, dx, boundary_distance
    logical :: local_ok, patch_ok, accepted, mpi_ok, physical_boundary
    integer :: nx, owner, ratio, subcycles, substep, child, child_count
    integer :: global_child
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
    patch_ok = .true.
    if (distribution%rank == owner) then
      allocate(state_start(solution%nvar, 0:nx + 1))
      allocate(state_end(solution%nvar, 0:nx + 1))
      allocate(flux(solution%nvar, 0:nx))
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
        left_integral = interval * flux(:, 0)
        right_integral = interval * flux(:, nx)
      end if
    end if
    call all_ranks_accept_sparse_1d( &
      distribution, patch_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) return
    call transfer_sparse_boundary_flux_integrals_1d( &
      distribution, owner, result_owner, left_integral, right_integral, &
      communication, local_ok)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) return
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
    local_ok = .true.
    if (distribution%rank == owner) then
      do child = 1, child_count
        call accumulate_coarse_flux_1d( &
          registers(level)%parents(parent_patch)%children(child), &
          flux(:, solution%hierarchy%relations(level)% &
            child_sets(parent_patch)%patches(child)%fine_coarse_lower - 1), &
          flux(:, solution%hierarchy%relations(level)% &
            child_sets(parent_patch)%patches(child)%fine_coarse_upper), &
          interval, local_ok)
        if (.not. local_ok) exit
      end do
    end if
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) return
    call distribute_sparse_interval_states_to_child_owners_1d( &
      distribution, solution, level, parent_patch, owner, state_start, &
      state_end, communication, local_ok)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) return
    allocate(child_left(solution%nvar, child_count))
    allocate(child_right(solution%nvar, child_count))
    child_left = 0.0_dp
    child_right = 0.0_dp

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
          child_right(:, child), owner, local_advances, communication, &
          local_ok)
        if (.not. local_ok) return
      end do
      call reconcile_sparse_adjacent_child_fluxes_1d( &
        species, distribution, solution, level, parent_patch, child_left, &
        child_right, communication, local_ok)
      call all_ranks_accept_sparse_1d( &
        distribution, local_ok, accepted, mpi_ok)
      if (.not. mpi_ok .or. .not. accepted) return
      local_ok = .true.
      if (distribution%rank == owner) then
        do child = 1, child_count
          call accumulate_fine_flux_1d( &
            registers(level)%parents(parent_patch)%children(child), &
            child_left(:, child) / child_interval, &
            child_right(:, child) / child_interval, child_interval, local_ok)
          if (.not. local_ok) exit
        end do
      end if
      call all_ranks_accept_sparse_1d( &
        distribution, local_ok, accepted, mpi_ok)
      if (.not. mpi_ok .or. .not. accepted) return
    end do

    call synchronize_sparse_parent_1d( &
      species, distribution, solution, registers, level, parent_patch, &
      local_ok)
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    ok = mpi_ok .and. accepted
  end subroutine advance_sparse_patch_transport_recursive_1d

  subroutine transfer_sparse_boundary_flux_integrals_1d( &
      distribution, source_owner, result_owner, left_integral, &
      right_integral, communication, ok)
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    integer, intent(in) :: source_owner, result_owner
    real(dp), intent(inout) :: left_integral(:), right_integral(:)
    type(mpi_amr_sparse_communication_counts_1d), intent(inout) :: &
      communication
    logical, intent(out) :: ok

    real(dp), allocatable :: payload(:)
    type(MPI_Status) :: status
    integer :: nvar, ierr

    ok = .false.
    nvar = size(left_integral)
    if (nvar < 1 .or. size(right_integral) /= nvar) return
    if (source_owner < 0 .or. source_owner >= distribution%nranks .or. &
        result_owner < 0 .or. result_owner >= distribution%nranks) return
    if (source_owner == result_owner) then
      ok = .true.
      return
    end if
    if (distribution%rank /= source_owner .and. &
        distribution%rank /= result_owner) then
      ok = .true.
      return
    end if

    allocate(payload(2 * nvar))
    if (distribution%rank == source_owner) then
      payload(1:nvar) = left_integral
      payload(nvar + 1:2 * nvar) = right_integral
      call MPI_Send( &
        payload, size(payload), MPI_DOUBLE_PRECISION, result_owner, &
        sparse_boundary_flux_tag, distribution%comm, ierr)
      if (ierr /= MPI_SUCCESS) return
      communication%boundary_flux_transfers = &
        communication%boundary_flux_transfers + 1
    else
      call MPI_Recv( &
        payload, size(payload), MPI_DOUBLE_PRECISION, source_owner, &
        sparse_boundary_flux_tag, distribution%comm, status, ierr)
      if (ierr /= MPI_SUCCESS) return
      left_integral = payload(1:nvar)
      right_integral = payload(nvar + 1:2 * nvar)
    end if
    ok = .true.
  end subroutine transfer_sparse_boundary_flux_integrals_1d

  subroutine distribute_sparse_interval_states_to_child_owners_1d( &
      distribution, solution, relation, parent, source_owner, state_start, &
      state_end, communication, ok)
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(in) :: solution
    integer, intent(in) :: relation, parent, source_owner
    real(dp), allocatable, intent(inout) :: state_start(:, :), state_end(:, :)
    type(mpi_amr_sparse_communication_counts_1d), intent(inout) :: &
      communication
    logical, intent(out) :: ok

    real(dp), allocatable :: payload(:)
    type(MPI_Status) :: status
    logical, allocatable :: recipients(:)
    logical :: has_remote
    integer :: child, child_count, child_index, child_owner
    integer :: nx, recipient, value_count, offset, ierr

    ok = .false.
    child_count = solution%hierarchy%relations(relation)% &
      child_sets(parent)%patch_count()
    nx = distribution%levels(relation)%cell_counts(parent)
    if (child_count < 1 .or. nx < 1) return
    allocate(recipients(distribution%nranks))
    recipients = .false.
    do child = 1, child_count
      child_index = solution%hierarchy%relations(relation)% &
        child_index(parent, child)
      child_owner = distribution%owner_of(relation, child_index)
      recipients(child_owner + 1) = .true.
    end do
    has_remote = .false.
    do recipient = 0, distribution%nranks - 1
      if (recipient /= source_owner .and. recipients(recipient + 1)) &
        has_remote = .true.
    end do
    value_count = 2 * solution%nvar * (nx + 2)

    if (distribution%rank == source_owner) then
      if (.not. allocated(state_start) .or. .not. allocated(state_end)) return
      if (size(state_start) /= solution%nvar * (nx + 2) .or. &
          size(state_end) /= solution%nvar * (nx + 2)) return
      if (has_remote) then
        allocate(payload(value_count))
        offset = 0
        call append_sparse_payload_2d(state_start, payload, offset)
        call append_sparse_payload_2d(state_end, payload, offset)
      end if
      do recipient = 0, distribution%nranks - 1
        if (recipient == source_owner .or. .not. recipients(recipient + 1)) &
          cycle
        call MPI_Send( &
          payload, value_count, MPI_DOUBLE_PRECISION, recipient, &
          sparse_interval_state_tag, distribution%comm, ierr)
        if (ierr /= MPI_SUCCESS) return
        communication%interval_state_transfers = &
          communication%interval_state_transfers + 1
      end do
    else if (recipients(distribution%rank + 1)) then
      allocate(state_start(solution%nvar, 0:nx + 1))
      allocate(state_end(solution%nvar, 0:nx + 1))
      allocate(payload(value_count))
      call MPI_Recv( &
        payload, value_count, MPI_DOUBLE_PRECISION, source_owner, &
        sparse_interval_state_tag, distribution%comm, status, ierr)
      if (ierr /= MPI_SUCCESS) return
      offset = 0
      call consume_sparse_payload_2d(payload, offset, state_start)
      call consume_sparse_payload_2d(payload, offset, state_end)
    end if
    ok = .true.
  end subroutine distribute_sparse_interval_states_to_child_owners_1d

  subroutine synchronize_sparse_counter_delta_1d( &
      distribution, baseline, counter, ok)
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    integer, intent(in) :: baseline(:)
    integer, intent(inout) :: counter(:)
    logical, intent(out) :: ok

    integer, allocatable :: local_delta(:), global_delta(:)
    integer :: ierr

    ok = .false.
    if (size(counter) /= size(baseline) .or. size(counter) < 1) return
    allocate(local_delta(size(counter)), global_delta(size(counter)))
    local_delta = counter - baseline
    if (any(local_delta < 0)) return
    call MPI_Allreduce( &
      local_delta, global_delta, size(counter), MPI_INTEGER, MPI_SUM, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    counter = baseline + global_delta
    ok = all(global_delta >= 0)
  end subroutine synchronize_sparse_counter_delta_1d

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
    integer :: child, child_count, child_index
    integer :: parent_owner, parent_nx

    ok = .false.
    child_count = solution%hierarchy%relations(level)% &
      child_sets(parent_patch)%patch_count()
    parent_owner = distribution%owner_of(level - 1, parent_patch)
    allocate(child_fields(child_count))
    do child = 1, child_count
      child_index = solution%hierarchy%relations(level)% &
        child_index(parent_patch, child)
      call transfer_sparse_child_interior_to_parent_1d( &
        distribution, solution, level, child_index, parent_owner, &
        child_state, local_ok)
      if (.not. local_ok) return
      if (distribution%rank == parent_owner) &
        child_fields(child)%values = child_state
      if (allocated(child_state)) deallocate(child_state)
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
      right_integrals, communication, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(inout) :: solution
    integer, intent(in) :: relation, parent
    real(dp), intent(inout) :: left_integrals(:, :), right_integrals(:, :)
    type(mpi_amr_sparse_communication_counts_1d), intent(inout) :: &
      communication
    logical, intent(out) :: ok

    real(dp), allocatable :: shared_integral(:), left_correction(:)
    real(dp), allocatable :: right_correction(:)
    logical, allocatable :: touched(:)
    logical :: local_ok
    real(dp) :: dx
    integer :: child, child_count, left_index, right_index, left_nx
    integer :: parent_owner

    ok = .false.
    child_count = solution%hierarchy%relations(relation)% &
      child_sets(parent)%patch_count()
    if (size(left_integrals, 1) /= solution%nvar .or. &
        size(right_integrals, 1) /= solution%nvar .or. &
        size(left_integrals, 2) /= child_count .or. &
        size(right_integrals, 2) /= child_count) return
    allocate(shared_integral(solution%nvar), left_correction(solution%nvar))
    allocate(right_correction(solution%nvar), touched(child_count))
    touched = .false.
    parent_owner = distribution%owner_of(relation - 1, parent)
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
      left_correction = 0.0_dp
      right_correction = 0.0_dp
      if (distribution%rank == parent_owner) then
        shared_integral = 0.5_dp * ( &
          right_integrals(:, child) + left_integrals(:, child + 1))
        left_correction = &
          -(shared_integral - right_integrals(:, child)) / dx
        right_correction = &
          (shared_integral - left_integrals(:, child + 1)) / dx
        right_integrals(:, child) = shared_integral
        left_integrals(:, child + 1) = shared_integral
      end if
      left_nx = distribution%levels(relation + 1)%cell_counts(left_index)
      call transfer_and_apply_sparse_flux_correction_1d( &
        distribution, solution, relation, left_index, left_nx, parent_owner, &
        left_correction, communication, local_ok)
      if (.not. local_ok) return
      call transfer_and_apply_sparse_flux_correction_1d( &
        distribution, solution, relation, right_index, 1, parent_owner, &
        right_correction, communication, local_ok)
      if (.not. local_ok) return
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

  subroutine transfer_and_apply_sparse_flux_correction_1d( &
      distribution, solution, relation, patch, cell, source_owner, &
      correction, communication, ok)
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(inout) :: solution
    integer, intent(in) :: relation, patch, cell, source_owner
    real(dp), intent(inout) :: correction(:)
    type(mpi_amr_sparse_communication_counts_1d), intent(inout) :: &
      communication
    logical, intent(out) :: ok

    type(MPI_Status) :: status
    integer :: patch_owner, ierr

    ok = .false.
    if (size(correction) /= solution%nvar) return
    patch_owner = distribution%owner_of(relation, patch)
    if (patch_owner == source_owner) then
      if (distribution%rank == patch_owner) &
        solution%levels(relation + 1)%patches(patch)%state(:, cell) = &
          solution%levels(relation + 1)%patches(patch)%state(:, cell) + &
          correction
      ok = .true.
      return
    end if

    if (distribution%rank == source_owner) then
      call MPI_Send( &
        correction, size(correction), MPI_DOUBLE_PRECISION, patch_owner, &
        sparse_flux_correction_tag, distribution%comm, ierr)
      if (ierr /= MPI_SUCCESS) return
      communication%shared_flux_correction_transfers = &
        communication%shared_flux_correction_transfers + 1
    else if (distribution%rank == patch_owner) then
      call MPI_Recv( &
        correction, size(correction), MPI_DOUBLE_PRECISION, source_owner, &
        sparse_flux_correction_tag, distribution%comm, status, ierr)
      if (ierr /= MPI_SUCCESS) return
      solution%levels(relation + 1)%patches(patch)%state(:, cell) = &
        solution%levels(relation + 1)%patches(patch)%state(:, cell) + &
        correction
    end if
    ok = .true.
  end subroutine transfer_and_apply_sparse_flux_correction_1d

  subroutine average_down_sparse_reactive_solution_1d( &
      species, distribution, solution, ok, local_parent_transfers)
    type(nasa7_species), intent(in) :: species(:)
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok
    integer, intent(inout), optional :: local_parent_transfers

    type(amr_level_field_1d), allocatable :: children(:)
    real(dp), allocatable :: child_state(:, :)
    logical :: local_ok, accepted, mpi_ok
    integer :: relation, parent, child, child_index, child_count
    integer :: parent_owner, parent_nx

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
          call transfer_sparse_child_interior_to_parent_1d( &
            distribution, solution, relation, child_index, parent_owner, &
            child_state, local_ok, local_parent_transfers)
          if (.not. local_ok) return
          if (distribution%rank == parent_owner) &
            children(child)%values = child_state
          if (allocated(child_state)) deallocate(child_state)
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

  subroutine transfer_sparse_child_interior_to_parent_1d( &
      distribution, solution, relation, child, parent_owner, child_state, &
      ok, local_parent_transfers)
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(in) :: solution
    integer, intent(in) :: relation, child, parent_owner
    real(dp), allocatable, intent(out) :: child_state(:, :)
    logical, intent(out) :: ok
    integer, intent(inout), optional :: local_parent_transfers

    type(MPI_Status) :: status
    integer :: child_owner, child_nx, ierr

    ok = .false.
    child_owner = distribution%owner_of(relation, child)
    child_nx = distribution%levels(relation + 1)%cell_counts(child)
    if (child_nx < 1) return
    if (child_owner == parent_owner) then
      if (distribution%rank == parent_owner) then
        allocate(child_state(solution%nvar, child_nx))
        child_state = solution%levels(relation + 1)% &
          patches(child)%state(:, 1:child_nx)
      end if
      ok = .true.
      return
    end if

    if (distribution%rank == child_owner) then
      call MPI_Send( &
        solution%levels(relation + 1)%patches(child)%state(:, 1:child_nx), &
        solution%nvar * child_nx, MPI_DOUBLE_PRECISION, parent_owner, &
        sparse_child_parent_tag, distribution%comm, ierr)
      if (ierr /= MPI_SUCCESS) return
      if (present(local_parent_transfers)) &
        local_parent_transfers = local_parent_transfers + 1
    else if (distribution%rank == parent_owner) then
      allocate(child_state(solution%nvar, child_nx))
      call MPI_Recv( &
        child_state, size(child_state), MPI_DOUBLE_PRECISION, child_owner, &
        sparse_child_parent_tag, distribution%comm, status, ierr)
      if (ierr /= MPI_SUCCESS) return
    end if
    ok = .true.
  end subroutine transfer_sparse_child_interior_to_parent_1d

  subroutine refresh_sparse_reactive_ghosts_1d( &
      species, config, distribution, solution, ok, local_halo_transfers, &
      local_parent_state_transfers)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok
    integer, intent(inout), optional :: local_halo_transfers
    integer, intent(inout), optional :: local_parent_state_transfers

    real(dp), allocatable :: parent_state(:, :)
    logical :: local_ok, accepted, mpi_ok
    integer :: relation, parent, child, child_index, child_count
    integer :: parent_owner, parent_nx

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
        call distribute_sparse_parent_state_to_child_owners_1d( &
          distribution, solution, relation, parent, parent_owner, &
          parent_state, local_ok, local_parent_state_transfers)
        if (.not. local_ok) return
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
        if (allocated(parent_state)) deallocate(parent_state)
        call all_ranks_accept_sparse_1d( &
          distribution, local_ok, accepted, mpi_ok)
        if (.not. mpi_ok .or. .not. accepted) return
        call exchange_sparse_adjacent_child_ghosts_1d( &
          config, distribution, solution, relation, parent, local_ok, &
          local_halo_transfers)
        call all_ranks_accept_sparse_1d( &
          distribution, local_ok, accepted, mpi_ok)
        if (.not. mpi_ok .or. .not. accepted) return
      end do
    end do
    ok = .true.
  end subroutine refresh_sparse_reactive_ghosts_1d

  subroutine distribute_sparse_parent_state_to_child_owners_1d( &
      distribution, solution, relation, parent, parent_owner, parent_state, &
      ok, local_parent_state_transfers)
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(in) :: solution
    integer, intent(in) :: relation, parent, parent_owner
    real(dp), allocatable, intent(out) :: parent_state(:, :)
    logical, intent(out) :: ok
    integer, intent(inout), optional :: local_parent_state_transfers

    type(MPI_Status) :: status
    logical, allocatable :: recipients(:)
    integer :: child, child_count, child_index, child_owner
    integer :: parent_nx, recipient, ierr

    ok = .false.
    child_count = solution%hierarchy%relations(relation)% &
      child_sets(parent)%patch_count()
    parent_nx = distribution%levels(relation)%cell_counts(parent)
    if (child_count < 1 .or. parent_nx < 1) return
    allocate(recipients(distribution%nranks))
    recipients = .false.
    do child = 1, child_count
      child_index = solution%hierarchy%relations(relation)% &
        child_index(parent, child)
      child_owner = distribution%owner_of(relation, child_index)
      recipients(child_owner + 1) = .true.
    end do

    if (distribution%rank == parent_owner) then
      allocate(parent_state(solution%nvar, 0:parent_nx + 1))
      parent_state = solution%levels(relation)%patches(parent)%state
      do recipient = 0, distribution%nranks - 1
        if (recipient == parent_owner .or. .not. recipients(recipient + 1)) &
          cycle
        call MPI_Send( &
          parent_state, size(parent_state), MPI_DOUBLE_PRECISION, recipient, &
          sparse_parent_state_tag, distribution%comm, ierr)
        if (ierr /= MPI_SUCCESS) return
        if (present(local_parent_state_transfers)) &
          local_parent_state_transfers = local_parent_state_transfers + 1
      end do
    else if (recipients(distribution%rank + 1)) then
      allocate(parent_state(solution%nvar, 0:parent_nx + 1))
      call MPI_Recv( &
        parent_state, size(parent_state), MPI_DOUBLE_PRECISION, parent_owner, &
        sparse_parent_state_tag, distribution%comm, status, ierr)
      if (ierr /= MPI_SUCCESS) return
    end if
    ok = .true.
  end subroutine distribute_sparse_parent_state_to_child_owners_1d

  subroutine exchange_sparse_adjacent_child_ghosts_1d( &
      config, distribution, solution, relation, parent, ok, &
      local_halo_transfers)
    type(reactive_1d_config), intent(in) :: config
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(inout) :: solution
    integer, intent(in) :: relation, parent
    logical, intent(out) :: ok
    integer, intent(inout), optional :: local_halo_transfers

    real(dp), allocatable :: send_payload(:), receive_payload(:)
    type(MPI_Status) :: status
    logical :: adjacent, wide_ghosts
    integer :: child, child_count, left_index, right_index
    integer :: left_owner, right_owner, left_nx, right_nx
    integer :: width, payload_size, layer, source_cell, offset, peer, ierr

    ok = .false.
    child_count = solution%hierarchy%relations(relation)% &
      child_sets(parent)%patch_count()
    wide_ghosts = sparse_uses_wide_ghosts(config)
    width = 1
    if (wide_ghosts) width = solution%ghost_width
    if (width < 1) return
    payload_size = (solution%nvar + 1) * width

    do child = 1, child_count - 1
      adjacent = solution%hierarchy%relations(relation)% &
        child_sets(parent)%patches(child)%fine%upper + 1 == &
        solution%hierarchy%relations(relation)% &
          child_sets(parent)%patches(child + 1)%fine%lower
      if (.not. adjacent) cycle
      left_index = solution%hierarchy%relations(relation)% &
        child_index(parent, child)
      right_index = solution%hierarchy%relations(relation)% &
        child_index(parent, child + 1)
      left_owner = distribution%owner_of(relation, left_index)
      right_owner = distribution%owner_of(relation, right_index)
      left_nx = distribution%levels(relation + 1)%cell_counts(left_index)
      right_nx = distribution%levels(relation + 1)%cell_counts(right_index)
      if (left_nx < width .or. right_nx < width) return

      if (left_owner == right_owner) then
        if (distribution%rank == left_owner) &
          call copy_local_sparse_adjacent_ghosts_1d( &
            solution, relation, left_index, right_index, left_nx, right_nx, &
            width, wide_ghosts)
        cycle
      end if
      if (distribution%rank /= left_owner .and. &
          distribution%rank /= right_owner) cycle

      allocate(send_payload(payload_size), receive_payload(payload_size))
      if (distribution%rank == left_owner) then
        peer = right_owner
        do layer = 1, width
          source_cell = left_nx - layer + 1
          offset = (layer - 1) * (solution%nvar + 1)
          send_payload(offset + 1:offset + solution%nvar) = &
            solution%levels(relation + 1)%patches(left_index)% &
              state(:, source_cell)
          send_payload(offset + solution%nvar + 1) = &
            solution%levels(relation + 1)%patches(left_index)% &
              temperature(source_cell)
        end do
      else
        peer = left_owner
        do layer = 1, width
          offset = (layer - 1) * (solution%nvar + 1)
          send_payload(offset + 1:offset + solution%nvar) = &
            solution%levels(relation + 1)%patches(right_index)%state(:, layer)
          send_payload(offset + solution%nvar + 1) = &
            solution%levels(relation + 1)%patches(right_index)% &
              temperature(layer)
        end do
      end if
      call MPI_Sendrecv( &
        send_payload, payload_size, MPI_DOUBLE_PRECISION, peer, &
        sparse_adjacent_halo_tag, receive_payload, payload_size, &
        MPI_DOUBLE_PRECISION, peer, sparse_adjacent_halo_tag, &
        distribution%comm, status, ierr)
      if (ierr /= MPI_SUCCESS) return

      if (distribution%rank == left_owner) then
        call unpack_sparse_right_ghosts_1d( &
          solution, relation, left_index, left_nx, width, wide_ghosts, &
          receive_payload)
        if (present(local_halo_transfers)) &
          local_halo_transfers = local_halo_transfers + 1
      else
        call unpack_sparse_left_ghosts_1d( &
          solution, relation, right_index, width, wide_ghosts, &
          receive_payload)
      end if
      deallocate(send_payload, receive_payload)
    end do
    ok = .true.
  end subroutine exchange_sparse_adjacent_child_ghosts_1d

  subroutine copy_local_sparse_adjacent_ghosts_1d( &
      solution, relation, left_index, right_index, left_nx, right_nx, &
      width, wide_ghosts)
    type(mpi_amr_sparse_reactive_solution_1d), intent(inout) :: solution
    integer, intent(in) :: relation, left_index, right_index
    integer, intent(in) :: left_nx, right_nx, width
    logical, intent(in) :: wide_ghosts

    integer :: layer

    solution%levels(relation + 1)%patches(right_index)%state(:, 0) = &
      solution%levels(relation + 1)%patches(left_index)%state(:, left_nx)
    solution%levels(relation + 1)%patches(right_index)%temperature(0) = &
      solution%levels(relation + 1)%patches(left_index)% &
        temperature(left_nx)
    solution%levels(relation + 1)%patches(left_index)% &
      state(:, left_nx + 1) = &
        solution%levels(relation + 1)%patches(right_index)%state(:, 1)
    solution%levels(relation + 1)%patches(left_index)% &
      temperature(left_nx + 1) = &
        solution%levels(relation + 1)%patches(right_index)%temperature(1)
    if (.not. wide_ghosts) return
    do layer = 1, width
      solution%levels(relation + 1)%patches(right_index)% &
        left_ghost_state(:, layer) = solution%levels(relation + 1)% &
          patches(left_index)%state(:, left_nx - layer + 1)
      solution%levels(relation + 1)%patches(right_index)% &
        left_ghost_temperature(layer) = solution%levels(relation + 1)% &
          patches(left_index)%temperature(left_nx - layer + 1)
      solution%levels(relation + 1)%patches(left_index)% &
        right_ghost_state(:, layer) = solution%levels(relation + 1)% &
          patches(right_index)%state(:, layer)
      solution%levels(relation + 1)%patches(left_index)% &
        right_ghost_temperature(layer) = solution%levels(relation + 1)% &
          patches(right_index)%temperature(layer)
    end do
  end subroutine copy_local_sparse_adjacent_ghosts_1d

  subroutine unpack_sparse_left_ghosts_1d( &
      solution, relation, patch, width, wide_ghosts, payload)
    type(mpi_amr_sparse_reactive_solution_1d), intent(inout) :: solution
    integer, intent(in) :: relation, patch, width
    logical, intent(in) :: wide_ghosts
    real(dp), intent(in) :: payload(:)

    integer :: layer, offset

    do layer = 1, width
      offset = (layer - 1) * (solution%nvar + 1)
      if (layer == 1) then
        solution%levels(relation + 1)%patches(patch)%state(:, 0) = &
          payload(offset + 1:offset + solution%nvar)
        solution%levels(relation + 1)%patches(patch)%temperature(0) = &
          payload(offset + solution%nvar + 1)
      end if
      if (wide_ghosts) then
        solution%levels(relation + 1)%patches(patch)% &
          left_ghost_state(:, layer) = &
            payload(offset + 1:offset + solution%nvar)
        solution%levels(relation + 1)%patches(patch)% &
          left_ghost_temperature(layer) = &
            payload(offset + solution%nvar + 1)
      end if
    end do
  end subroutine unpack_sparse_left_ghosts_1d

  subroutine unpack_sparse_right_ghosts_1d( &
      solution, relation, patch, nx, width, wide_ghosts, payload)
    type(mpi_amr_sparse_reactive_solution_1d), intent(inout) :: solution
    integer, intent(in) :: relation, patch, nx, width
    logical, intent(in) :: wide_ghosts
    real(dp), intent(in) :: payload(:)

    integer :: layer, offset

    do layer = 1, width
      offset = (layer - 1) * (solution%nvar + 1)
      if (layer == 1) then
        solution%levels(relation + 1)%patches(patch)%state(:, nx + 1) = &
          payload(offset + 1:offset + solution%nvar)
        solution%levels(relation + 1)%patches(patch)% &
          temperature(nx + 1) = payload(offset + solution%nvar + 1)
      end if
      if (wide_ghosts) then
        solution%levels(relation + 1)%patches(patch)% &
          right_ghost_state(:, layer) = &
            payload(offset + 1:offset + solution%nvar)
        solution%levels(relation + 1)%patches(patch)% &
          right_ghost_temperature(layer) = &
            payload(offset + solution%nvar + 1)
      end if
    end do
  end subroutine unpack_sparse_right_ghosts_1d

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
      old_distribution, new_distribution, old_sparse, new_sparse, ok, &
      local_patch_transfers)
    type(mpi_amr_patch_distribution_1d), intent(in) :: old_distribution
    type(mpi_amr_patch_distribution_1d), intent(in) :: new_distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(in) :: old_sparse
    type(mpi_amr_sparse_reactive_solution_1d), intent(out) :: new_sparse
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_patch_transfers

    logical :: local_ok, accepted, mpi_ok
    integer :: level, patch, nx, old_owner, new_owner, transfers

    ok = .false.
    transfers = 0
    if (present(local_patch_transfers)) local_patch_transfers = 0
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
        call all_ranks_accept_sparse_1d( &
          old_distribution, local_ok, accepted, mpi_ok)
        if (.not. mpi_ok .or. .not. accepted) then
          new_sparse = mpi_amr_sparse_reactive_solution_1d()
          return
        end if
        if (old_owner /= new_owner .and. &
            old_distribution%rank == old_owner) transfers = transfers + 1
      end do
    end do
    local_ok = new_sparse%is_valid(new_distribution)
    call all_ranks_accept_sparse_1d( &
      old_distribution, local_ok, accepted, mpi_ok)
    ok = mpi_ok .and. accepted
    if (.not. ok) new_sparse = mpi_amr_sparse_reactive_solution_1d()
    if (ok .and. present(local_patch_transfers)) &
      local_patch_transfers = transfers
  end subroutine migrate_owned_patch_tree_reactive_1d

  subroutine regrid_sparse_patch_tree_reactive_1d( &
      species, config, plans, old_distribution, solution, new_distribution, &
      changed, transferred_cells, ok, local_prolongation_transfers, &
      local_overlap_transfers)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    type(amr_patch_level_plan_1d), intent(in) :: plans(:)
    type(mpi_amr_patch_distribution_1d), intent(in) :: old_distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(inout) :: solution
    type(mpi_amr_patch_distribution_1d), intent(out) :: new_distribution
    logical, intent(out) :: changed
    integer, intent(out) :: transferred_cells
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_prolongation_transfers
    integer, intent(out), optional :: local_overlap_transfers

    type(mpi_amr_sparse_reactive_solution_1d) :: backup
    type(mpi_amr_sparse_reactive_solution_1d) :: rebuilt
    type(mpi_amr_patch_distribution_1d) :: rebuilt_distribution
    type(amr_patch_tree_hierarchy_1d) :: rebuilt_hierarchy
    real(dp) :: tolerance
    logical :: local_ok, accepted, mpi_ok
    integer :: common_levels, prolongation_transfers, overlap_transfers

    ok = .false.
    changed = .false.
    transferred_cells = 0
    prolongation_transfers = 0
    overlap_transfers = 0
    if (present(local_prolongation_transfers)) &
      local_prolongation_transfers = 0
    if (present(local_overlap_transfers)) local_overlap_transfers = 0
    new_distribution = old_distribution
    local_ok = size(species) >= 1 .and. config%amr_enabled .and. &
      solution%is_valid(old_distribution) .and. &
      config%nx == solution%hierarchy%base_cells
    if (local_ok) then
      tolerance = 128.0_dp * epsilon(1.0_dp) * max( &
        1.0_dp, abs(config%x_lower), abs(config%x_upper), &
        abs(solution%hierarchy%x_lower), abs(solution%hierarchy%x_upper))
      local_ok = abs(config%x_lower - solution%hierarchy%x_lower) <= &
          tolerance .and. &
        abs(config%x_upper - solution%hierarchy%x_upper) <= tolerance
    end if
    call all_ranks_accept_sparse_1d( &
      old_distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) return
    backup = solution

    call initialize_patch_tree_1d( &
      config%nx, config%x_lower, config%x_upper, plans, &
      rebuilt_hierarchy, local_ok)
    if (local_ok) local_ok = &
      sparse_patch_tree_children_are_interior_1d(rebuilt_hierarchy)
    call all_ranks_accept_sparse_1d( &
      old_distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) go to 900

    changed = .not. same_sparse_patch_tree_hierarchy_1d( &
      solution%hierarchy, rebuilt_hierarchy)
    if (.not. changed) then
      solution%regrid_evaluations = solution%regrid_evaluations + 1
      local_ok = solution%is_valid(old_distribution)
      call all_ranks_accept_sparse_1d( &
        old_distribution, local_ok, accepted, mpi_ok)
      if (.not. mpi_ok .or. .not. accepted) go to 900
      ok = .true.
      return
    end if

    call average_down_sparse_reactive_solution_1d( &
      species, old_distribution, solution, local_ok)
    call all_ranks_accept_sparse_1d( &
      old_distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) go to 900

    call initialize_mpi_amr_patch_distribution_1d( &
      rebuilt_hierarchy, old_distribution%comm, rebuilt_distribution, &
      local_ok)
    call all_ranks_accept_sparse_1d( &
      old_distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) go to 900
    call initialize_direct_sparse_regrid_1d( &
      species, old_distribution, rebuilt_distribution, solution, &
      rebuilt_hierarchy, rebuilt, local_ok, prolongation_transfers)
    call all_ranks_accept_sparse_1d( &
      old_distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) go to 900
    call transfer_direct_sparse_regrid_overlap_1d( &
      old_distribution, rebuilt_distribution, solution, rebuilt, &
      transferred_cells, local_ok, overlap_transfers)
    call all_ranks_accept_sparse_1d( &
      old_distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) go to 900
    call average_down_sparse_reactive_solution_1d( &
      species, rebuilt_distribution, rebuilt, local_ok)
    call all_ranks_accept_sparse_1d( &
      old_distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) go to 900
    call refresh_sparse_reactive_ghosts_1d( &
      species, config, rebuilt_distribution, rebuilt, local_ok)
    call all_ranks_accept_sparse_1d( &
      old_distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) go to 900

    rebuilt%time = backup%time
    rebuilt%steps = backup%steps
    common_levels = min(size(backup%levels), size(rebuilt%levels))
    rebuilt%level_advances(1:common_levels) = &
      backup%level_advances(1:common_levels)
    rebuilt%transport_level_advances(1:common_levels) = &
      backup%transport_level_advances(1:common_levels)
    rebuilt%regrid_evaluations = backup%regrid_evaluations + 1
    rebuilt%regrids = backup%regrids + 1
    rebuilt%overlap_cells_transferred = &
      backup%overlap_cells_transferred + transferred_cells
    local_ok = rebuilt%is_valid(rebuilt_distribution)
    call all_ranks_accept_sparse_1d( &
      old_distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) go to 900
    solution = rebuilt
    new_distribution = rebuilt_distribution
    if (present(local_prolongation_transfers)) &
      local_prolongation_transfers = prolongation_transfers
    if (present(local_overlap_transfers)) &
      local_overlap_transfers = overlap_transfers
    ok = .true.
    return

900 continue
    solution = backup
    new_distribution = old_distribution
    changed = .false.
    transferred_cells = 0
    if (present(local_prolongation_transfers)) &
      local_prolongation_transfers = 0
    if (present(local_overlap_transfers)) local_overlap_transfers = 0
    ok = .false.
  end subroutine regrid_sparse_patch_tree_reactive_1d

  subroutine initialize_direct_sparse_regrid_1d( &
      species, old_distribution, new_distribution, old_solution, hierarchy, &
      rebuilt, ok, local_prolongation_transfers)
    type(nasa7_species), intent(in) :: species(:)
    type(mpi_amr_patch_distribution_1d), intent(in) :: old_distribution
    type(mpi_amr_patch_distribution_1d), intent(in) :: new_distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(in) :: old_solution
    type(amr_patch_tree_hierarchy_1d), intent(in) :: hierarchy
    type(mpi_amr_sparse_reactive_solution_1d), intent(out) :: rebuilt
    logical, intent(out) :: ok
    integer, intent(out) :: local_prolongation_transfers

    logical :: local_ok, accepted, mpi_ok
    integer :: level, relation, parent, patch_count
    integer :: old_root_owner, new_root_owner, root_nx

    ok = .false.
    local_prolongation_transfers = 0
    rebuilt%hierarchy = hierarchy
    rebuilt%rank = new_distribution%rank
    rebuilt%nranks = new_distribution%nranks
    rebuilt%nvar = old_solution%nvar
    rebuilt%ghost_width = old_solution%ghost_width
    allocate(rebuilt%level_advances(hierarchy%level_count()))
    allocate(rebuilt%transport_level_advances(hierarchy%level_count()))
    rebuilt%level_advances = 0
    rebuilt%transport_level_advances = 0
    allocate(rebuilt%levels(hierarchy%level_count()))
    do level = 1, size(rebuilt%levels)
      patch_count = hierarchy%level_patch_count(level - 1)
      allocate(rebuilt%levels(level)%patches(patch_count))
      allocate(rebuilt%levels(level)%is_local(patch_count))
      rebuilt%levels(level)%is_local = &
        new_distribution%levels(level)%owners == new_distribution%rank
    end do

    old_root_owner = old_distribution%owner_of(0, 1)
    new_root_owner = new_distribution%owner_of(0, 1)
    root_nx = new_distribution%levels(1)%cell_counts(1)
    call migrate_one_patch_1d( &
      old_distribution%comm, old_distribution%rank, old_root_owner, &
      new_root_owner, rebuilt%nvar, root_nx, rebuilt%ghost_width, &
      old_solution%levels(1)%patches(1), &
      rebuilt%levels(1)%patches(1), local_ok)
    call all_ranks_accept_sparse_1d( &
      old_distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) return

    do relation = 1, size(hierarchy%relations)
      do parent = 1, hierarchy%relations(relation)%parent_patch_count()
        call prolong_direct_sparse_regrid_parent_1d( &
          species, new_distribution, rebuilt, relation, parent, local_ok, &
          local_prolongation_transfers)
        if (.not. local_ok) return
      end do
    end do
    ok = rebuilt%is_valid(new_distribution)
  end subroutine initialize_direct_sparse_regrid_1d

  subroutine prolong_direct_sparse_regrid_parent_1d( &
      species, distribution, solution, relation, parent, ok, &
      local_prolongation_transfers)
    type(nasa7_species), intent(in) :: species(:)
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(inout) :: solution
    integer, intent(in) :: relation, parent
    logical, intent(out) :: ok
    integer, intent(inout) :: local_prolongation_transfers

    type(amr_level_field_1d), allocatable :: children(:)
    type(MPI_Status) :: status
    logical :: local_ok, accepted, mpi_ok
    integer :: parent_owner, parent_nx, child, child_index, child_owner
    integer :: child_nx, ierr

    ok = .false.
    parent_owner = distribution%owner_of(relation - 1, parent)
    parent_nx = distribution%levels(relation)%cell_counts(parent)
    local_ok = parent_owner >= 0 .and. parent_nx >= 1
    if (local_ok .and. distribution%rank == parent_owner) then
      call prolong_patch_set_1d( &
        solution%levels(relation)%patches(parent)%state(:, 1:parent_nx), &
        solution%hierarchy%relations(relation)%child_sets(parent), &
        children, local_ok)
    end if
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) return

    do child = 1, solution%hierarchy%relations(relation)% &
        child_sets(parent)%patch_count()
      child_index = solution%hierarchy%relations(relation)% &
        child_index(parent, child)
      child_owner = distribution%owner_of(relation, child_index)
      child_nx = distribution%levels(relation + 1)%cell_counts(child_index)
      local_ok = child_owner >= 0 .and. child_nx >= 1
      if (.not. local_ok) then
        call all_ranks_accept_sparse_1d( &
          distribution, local_ok, accepted, mpi_ok)
        return
      end if

      if (parent_owner == child_owner) then
        if (distribution%rank == child_owner) then
          call allocate_sparse_patch( &
            solution%levels(relation + 1)%patches(child_index), &
            solution%nvar, child_nx, solution%ghost_width)
          solution%levels(relation + 1)%patches(child_index)% &
            state(:, 1:child_nx) = children(child)%values
        end if
      else if (distribution%rank == parent_owner) then
        call MPI_Send( &
          children(child)%values, solution%nvar * child_nx, &
          MPI_DOUBLE_PRECISION, child_owner, sparse_regrid_prolongation_tag, &
          distribution%comm, ierr)
        local_ok = ierr == MPI_SUCCESS
        if (local_ok) local_prolongation_transfers = &
          local_prolongation_transfers + 1
      else if (distribution%rank == child_owner) then
        call allocate_sparse_patch( &
          solution%levels(relation + 1)%patches(child_index), &
          solution%nvar, child_nx, solution%ghost_width)
        call MPI_Recv( &
          solution%levels(relation + 1)%patches(child_index)% &
            state(:, 1:child_nx), &
          solution%nvar * child_nx, MPI_DOUBLE_PRECISION, parent_owner, &
          sparse_regrid_prolongation_tag, distribution%comm, status, ierr)
        local_ok = ierr == MPI_SUCCESS
      end if
      if (local_ok .and. distribution%rank == child_owner) &
        call recover_level_temperatures_1d( &
          species, solution%levels(relation + 1)%patches(child_index)%state, &
          solution%levels(relation + 1)%patches(child_index)%temperature, &
          child_nx, local_ok)
      call all_ranks_accept_sparse_1d( &
        distribution, local_ok, accepted, mpi_ok)
      if (.not. mpi_ok .or. .not. accepted) return
    end do
    ok = .true.
  end subroutine prolong_direct_sparse_regrid_parent_1d

  subroutine transfer_direct_sparse_regrid_overlap_1d( &
      old_distribution, new_distribution, old_solution, new_solution, &
      transferred_cells, ok, local_overlap_transfers)
    type(mpi_amr_patch_distribution_1d), intent(in) :: old_distribution
    type(mpi_amr_patch_distribution_1d), intent(in) :: new_distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(in) :: old_solution
    type(mpi_amr_sparse_reactive_solution_1d), intent(inout) :: new_solution
    integer, intent(out) :: transferred_cells
    logical, intent(out) :: ok
    integer, intent(out) :: local_overlap_transfers

    type(amr_two_level_hierarchy_1d) :: old_geometry, new_geometry
    real(dp) :: old_lower, old_upper, new_lower, new_upper
    real(dp) :: overlap_lower, overlap_upper, old_dx, new_dx, tolerance
    real(dp) :: old_offset, new_offset
    logical :: local_ok, accepted, mpi_ok
    integer :: level, old_patch, new_patch, common_levels
    integer :: old_first, new_first, cell_count
    integer :: old_owner, new_owner

    ok = .false.
    transferred_cells = 0
    local_overlap_transfers = 0
    common_levels = min(size(old_solution%levels), size(new_solution%levels))
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
        if (.not. local_ok) return
        call sparse_patch_physical_bounds_1d( &
          old_geometry, old_lower, old_upper)
        old_owner = old_distribution%owner_of(level - 1, old_patch)
        do new_patch = 1, size(new_solution%levels(level)%patches)
          call patch_tree_child_geometry_1d( &
            new_solution%hierarchy%relations(level - 1), new_patch, &
            new_geometry, local_ok)
          if (.not. local_ok) return
          call sparse_patch_physical_bounds_1d( &
            new_geometry, new_lower, new_upper)
          overlap_lower = max(old_lower, new_lower)
          overlap_upper = min(old_upper, new_upper)
          if (overlap_upper <= overlap_lower + tolerance) cycle
          old_offset = (overlap_lower - old_lower) / old_dx
          new_offset = (overlap_lower - new_lower) / new_dx
          local_ok = &
            abs(old_offset - real(nint(old_offset), dp)) <= &
              tolerance / old_dx .and. &
            abs(new_offset - real(nint(new_offset), dp)) <= &
              tolerance / new_dx
          if (.not. local_ok) return
          old_first = nint(old_offset) + 1
          new_first = nint(new_offset) + 1
          cell_count = nint((overlap_upper - overlap_lower) / old_dx)
          if (cell_count < 1) cycle
          new_owner = new_distribution%owner_of(level - 1, new_patch)
          call transfer_sparse_regrid_overlap_segment_1d( &
            old_distribution, old_owner, new_owner, old_solution%nvar, &
            old_solution%levels(level)%patches(old_patch), old_first, &
            new_solution%levels(level)%patches(new_patch), new_first, &
            cell_count, local_ok)
          call all_ranks_accept_sparse_1d( &
            old_distribution, local_ok, accepted, mpi_ok)
          if (.not. mpi_ok .or. .not. accepted) return
          transferred_cells = transferred_cells + cell_count
          if (old_owner /= new_owner .and. &
              old_distribution%rank == old_owner) &
            local_overlap_transfers = local_overlap_transfers + 1
        end do
      end do
    end do
    ok = .true.
  end subroutine transfer_direct_sparse_regrid_overlap_1d

  subroutine transfer_sparse_regrid_overlap_segment_1d( &
      distribution, old_owner, new_owner, nvar, old_patch, old_first, &
      new_patch, new_first, cell_count, ok)
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    integer, intent(in) :: old_owner, new_owner, nvar
    type(amr_patch_tree_reactive_patch_1d), intent(in) :: old_patch
    integer, intent(in) :: old_first
    type(amr_patch_tree_reactive_patch_1d), intent(inout) :: new_patch
    integer, intent(in) :: new_first, cell_count
    logical, intent(out) :: ok

    real(dp), allocatable :: payload(:)
    type(MPI_Status) :: status
    integer :: state_count, ierr

    ok = old_owner >= 0 .and. new_owner >= 0 .and. nvar >= 1 .and. &
      old_first >= 1 .and. new_first >= 1 .and. cell_count >= 1
    if (.not. ok) return
    if (old_owner == new_owner) then
      if (distribution%rank == old_owner) then
        new_patch%state(:, new_first:new_first + cell_count - 1) = &
          old_patch%state(:, old_first:old_first + cell_count - 1)
        new_patch%temperature(new_first:new_first + cell_count - 1) = &
          old_patch%temperature(old_first:old_first + cell_count - 1)
      end if
      return
    end if
    if (distribution%rank /= old_owner .and. &
        distribution%rank /= new_owner) return

    state_count = nvar * cell_count
    allocate(payload(state_count + cell_count))
    if (distribution%rank == old_owner) then
      payload(1:state_count) = reshape( &
        old_patch%state(:, old_first:old_first + cell_count - 1), &
        [state_count])
      payload(state_count + 1:) = &
        old_patch%temperature(old_first:old_first + cell_count - 1)
      call MPI_Send( &
        payload, size(payload), MPI_DOUBLE_PRECISION, new_owner, &
        sparse_regrid_overlap_tag, distribution%comm, ierr)
      ok = ierr == MPI_SUCCESS
      return
    end if

    call MPI_Recv( &
      payload, size(payload), MPI_DOUBLE_PRECISION, old_owner, &
      sparse_regrid_overlap_tag, distribution%comm, status, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    new_patch%state(:, new_first:new_first + cell_count - 1) = &
      reshape(payload(1:state_count), [nvar, cell_count])
    new_patch%temperature(new_first:new_first + cell_count - 1) = &
      payload(state_count + 1:)
    ok = .true.
  end subroutine transfer_sparse_regrid_overlap_segment_1d

  pure logical function same_sparse_patch_tree_hierarchy_1d(first, second) &
      result(same)
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
  end function same_sparse_patch_tree_hierarchy_1d

  pure logical function sparse_patch_tree_children_are_interior_1d( &
      hierarchy) result(interior)
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
  end function sparse_patch_tree_children_are_interior_1d

  pure subroutine sparse_patch_physical_bounds_1d(geometry, lower, upper)
    type(amr_two_level_hierarchy_1d), intent(in) :: geometry
    real(dp), intent(out) :: lower, upper

    lower = geometry%x_lower + &
      real(geometry%fine_coarse_lower - 1, dp) * geometry%coarse_dx
    upper = geometry%x_lower + &
      real(geometry%fine_coarse_upper, dp) * geometry%coarse_dx
  end subroutine sparse_patch_physical_bounds_1d

  subroutine regrid_tagged_sparse_patch_tree_reactive_1d( &
      species, config, old_distribution, solution, new_distribution, &
      changed, tagged_cells, transferred_cells, ok, &
      local_tagging_evaluations, local_candidate_transfers, &
      local_prolongation_transfers, local_overlap_transfers)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    type(mpi_amr_patch_distribution_1d), intent(in) :: old_distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(inout) :: solution
    type(mpi_amr_patch_distribution_1d), intent(out) :: new_distribution
    logical, intent(out) :: changed
    integer, intent(out) :: tagged_cells, transferred_cells
    logical, intent(out) :: ok
    integer, intent(out), optional :: local_tagging_evaluations
    integer, intent(out), optional :: local_candidate_transfers
    integer, intent(out), optional :: local_prolongation_transfers
    integer, intent(out), optional :: local_overlap_transfers

    type(mpi_amr_sparse_reactive_solution_1d) :: backup, planning_solution
    type(amr_patch_level_plan_1d), allocatable :: plans(:)
    logical :: local_ok, accepted, mpi_ok
    integer :: tagging_evaluations, candidate_transfers
    integer :: prolongation_transfers, overlap_transfers

    ok = .false.
    changed = .false.
    tagged_cells = 0
    transferred_cells = 0
    tagging_evaluations = 0
    candidate_transfers = 0
    prolongation_transfers = 0
    overlap_transfers = 0
    if (present(local_tagging_evaluations)) local_tagging_evaluations = 0
    if (present(local_candidate_transfers)) local_candidate_transfers = 0
    if (present(local_prolongation_transfers)) &
      local_prolongation_transfers = 0
    if (present(local_overlap_transfers)) local_overlap_transfers = 0
    new_distribution = old_distribution
    local_ok = size(species) >= 1 .and. &
      solution%is_valid(old_distribution)
    call all_ranks_accept_sparse_1d( &
      old_distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) return
    backup = solution
    planning_solution = solution

    call average_down_sparse_reactive_solution_1d( &
      species, old_distribution, planning_solution, local_ok)
    call all_ranks_accept_sparse_1d( &
      old_distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) go to 900
    call plan_tagged_sparse_patch_tree_reactive_1d( &
      species, config, old_distribution, planning_solution, plans, &
      tagged_cells, local_ok, tagging_evaluations, candidate_transfers)
    call all_ranks_accept_sparse_1d( &
      old_distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) go to 900
    call regrid_sparse_patch_tree_reactive_1d( &
      species, config, plans, old_distribution, solution, new_distribution, &
      changed, transferred_cells, local_ok, prolongation_transfers, &
      overlap_transfers)
    call all_ranks_accept_sparse_1d( &
      old_distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) go to 900
    if (present(local_tagging_evaluations)) &
      local_tagging_evaluations = tagging_evaluations
    if (present(local_candidate_transfers)) &
      local_candidate_transfers = candidate_transfers
    if (present(local_prolongation_transfers)) &
      local_prolongation_transfers = prolongation_transfers
    if (present(local_overlap_transfers)) &
      local_overlap_transfers = overlap_transfers
    ok = .true.
    return

900 continue
    solution = backup
    new_distribution = old_distribution
    changed = .false.
    tagged_cells = 0
    transferred_cells = 0
    if (present(local_tagging_evaluations)) local_tagging_evaluations = 0
    if (present(local_candidate_transfers)) local_candidate_transfers = 0
    if (present(local_prolongation_transfers)) &
      local_prolongation_transfers = 0
    if (present(local_overlap_transfers)) local_overlap_transfers = 0
    ok = .false.
  end subroutine regrid_tagged_sparse_patch_tree_reactive_1d

  subroutine plan_tagged_sparse_patch_tree_reactive_1d( &
      species, config, distribution, solution, plans, tagged_cells, ok, &
      local_tagging_evaluations, local_candidate_transfers)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    type(mpi_amr_patch_distribution_1d), intent(in) :: distribution
    type(mpi_amr_sparse_reactive_solution_1d), intent(in) :: solution
    type(amr_patch_level_plan_1d), allocatable, intent(out) :: plans(:)
    integer, intent(out) :: tagged_cells
    logical, intent(out) :: ok
    integer, intent(out) :: local_tagging_evaluations
    integer, intent(out) :: local_candidate_transfers

    type(amr_patch_level_plan_1d), allocatable :: workspace(:)
    type(amr_child_patch_plan_1d), allocatable :: extended_patches(:)
    type(amr_patch_tree_hierarchy_1d) :: candidate_hierarchy
    type(mpi_amr_patch_distribution_1d) :: candidate_distribution
    type(mpi_amr_sparse_reactive_solution_1d) :: candidate_solution
    type(amr_regrid_plan_collection_1d) :: collection
    integer, allocatable :: local_bounds(:), global_bounds(:)
    integer :: local_header(2), global_header(2)
    logical :: local_ok, accepted, mpi_ok
    integer :: maximum_relations, relation_count, relation
    integer :: parent_count, parent, parent_owner, parent_nx
    integer :: child_count, child, entry, offset, ierr, step_transfers

    ok = .false.
    tagged_cells = 0
    local_tagging_evaluations = 0
    local_candidate_transfers = 0
    local_ok = size(species) >= 1 .and. config%amr_max_levels >= 2 .and. &
      solution%is_valid(distribution) .and. config%nx == &
        solution%hierarchy%base_cells
    call all_ranks_accept_sparse_1d( &
      distribution, local_ok, accepted, mpi_ok)
    if (.not. mpi_ok .or. .not. accepted) return

    maximum_relations = config%amr_max_levels - 1
    allocate(workspace(maximum_relations))
    candidate_distribution = distribution
    candidate_solution = solution
    relation_count = 0
    do relation = 1, maximum_relations
      parent_count = candidate_solution%hierarchy% &
        level_patch_count(relation - 1)
      child_count = 0
      do parent = 1, parent_count
        parent_owner = candidate_distribution%owner_of(relation - 1, parent)
        parent_nx = candidate_distribution%levels(relation)% &
          cell_counts(parent)
        local_header = 0
        local_ok = parent_owner >= 0 .and. parent_nx >= 1
        if (local_ok .and. distribution%rank == parent_owner) then
          call plan_tagged_reactive_parent_1d( &
            config, candidate_solution%levels(relation)%patches(parent)% &
              state(:, 1:parent_nx), collection, local_ok)
          if (local_ok) then
            local_header = [collection%tagged_cell_count, &
              collection%patch_count()]
            local_tagging_evaluations = local_tagging_evaluations + 1
          end if
        end if
        call all_ranks_accept_sparse_1d( &
          distribution, local_ok, accepted, mpi_ok)
        if (.not. mpi_ok .or. .not. accepted) return
        call MPI_Allreduce( &
          local_header, global_header, 2, MPI_INTEGER, MPI_SUM, &
          distribution%comm, ierr)
        if (ierr /= MPI_SUCCESS .or. any(global_header < 0)) return
        tagged_cells = tagged_cells + global_header(1)
        if (global_header(2) == 0) cycle

        allocate(local_bounds(2 * global_header(2)))
        allocate(global_bounds(2 * global_header(2)))
        local_bounds = 0
        if (distribution%rank == parent_owner) then
          do child = 1, global_header(2)
            local_bounds(2 * child - 1) = &
              collection%plans(child)%patch_lower
            local_bounds(2 * child) = collection%plans(child)%patch_upper
          end do
        end if
        call MPI_Allreduce( &
          local_bounds, global_bounds, size(local_bounds), MPI_INTEGER, &
          MPI_SUM, distribution%comm, ierr)
        if (ierr /= MPI_SUCCESS) return
        offset = child_count
        child_count = child_count + global_header(2)
        allocate(extended_patches(child_count))
        if (offset > 0) extended_patches(1:offset) = &
          workspace(relation)%patches
        do child = 1, global_header(2)
          entry = offset + child
          extended_patches(entry)%parent_patch = parent
          extended_patches(entry)%lower = &
            global_bounds(2 * child - 1)
          extended_patches(entry)%upper = &
            global_bounds(2 * child)
        end do
        call move_alloc(extended_patches, workspace(relation)%patches)
        deallocate(local_bounds, global_bounds)
      end do
      if (child_count == 0) exit

      workspace(relation)%refinement_ratio = config%amr_refinement_ratio
      relation_count = relation
      call initialize_patch_tree_1d( &
        config%nx, config%x_lower, config%x_upper, &
        workspace(1:relation_count), candidate_hierarchy, local_ok)
      if (local_ok) local_ok = &
        sparse_patch_tree_children_are_interior_1d(candidate_hierarchy)
      call all_ranks_accept_sparse_1d( &
        distribution, local_ok, accepted, mpi_ok)
      if (.not. mpi_ok .or. .not. accepted) return
      call initialize_mpi_amr_patch_distribution_1d( &
        candidate_hierarchy, distribution%comm, candidate_distribution, &
        local_ok)
      call all_ranks_accept_sparse_1d( &
        distribution, local_ok, accepted, mpi_ok)
      if (.not. mpi_ok .or. .not. accepted) return
      call initialize_direct_sparse_regrid_1d( &
        species, distribution, candidate_distribution, solution, &
        candidate_hierarchy, candidate_solution, local_ok, step_transfers)
      call all_ranks_accept_sparse_1d( &
        distribution, local_ok, accepted, mpi_ok)
      if (.not. mpi_ok .or. .not. accepted) return
      local_candidate_transfers = &
        local_candidate_transfers + step_transfers
    end do

    allocate(plans(relation_count))
    if (relation_count > 0) plans = workspace(1:relation_count)
    ok = .true.
  end subroutine plan_tagged_sparse_patch_tree_reactive_1d

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

    real(dp), allocatable :: payload(:)
    type(MPI_Status) :: status
    integer :: value_count
    integer :: ierr

    ok = .true.
    if (old_owner < 0 .or. new_owner < 0 .or. nvar < 1 .or. nx < 1 .or. &
        ghost_width < 1) then
      ok = .false.
      return
    end if
    if (old_owner == new_owner) then
      if (rank == old_owner) destination = source
      return
    end if
    if (rank /= old_owner .and. rank /= new_owner) return

    value_count = nvar * (nx + 2) + (nx + 2) + &
      2 * nvar * ghost_width + 2 * ghost_width
    allocate(payload(value_count))
    if (rank == old_owner) then
      call pack_sparse_patch_values(source, payload, ok)
      if (.not. ok) return
      call MPI_Send(payload, value_count, MPI_DOUBLE_PRECISION, new_owner, &
        sparse_patch_migration_tag, comm, ierr)
      ok = ierr == MPI_SUCCESS
      return
    end if

    call MPI_Recv(payload, value_count, MPI_DOUBLE_PRECISION, old_owner, &
      sparse_patch_migration_tag, comm, status, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call allocate_sparse_patch(destination, nvar, nx, ghost_width)
    call unpack_sparse_patch_values(payload, destination, ok)
  end subroutine migrate_one_patch_1d

  subroutine pack_sparse_patch_values(patch, payload, ok)
    type(amr_patch_tree_reactive_patch_1d), intent(in) :: patch
    real(dp), intent(out) :: payload(:)
    logical, intent(out) :: ok

    integer :: offset, count

    ok = allocated(patch%state) .and. allocated(patch%temperature) .and. &
      allocated(patch%left_ghost_state) .and. &
      allocated(patch%right_ghost_state) .and. &
      allocated(patch%left_ghost_temperature) .and. &
      allocated(patch%right_ghost_temperature)
    if (.not. ok) return
    count = sparse_patch_value_count(patch)
    ok = size(payload) == count
    if (.not. ok) return
    offset = 0
    call append_sparse_payload_2d(patch%state, payload, offset)
    call append_sparse_payload_1d(patch%temperature, payload, offset)
    call append_sparse_payload_2d( &
      patch%left_ghost_state, payload, offset)
    call append_sparse_payload_2d( &
      patch%right_ghost_state, payload, offset)
    call append_sparse_payload_1d( &
      patch%left_ghost_temperature, payload, offset)
    call append_sparse_payload_1d( &
      patch%right_ghost_temperature, payload, offset)
    ok = offset == size(payload)
  end subroutine pack_sparse_patch_values

  subroutine unpack_sparse_patch_values(payload, patch, ok)
    real(dp), intent(in) :: payload(:)
    type(amr_patch_tree_reactive_patch_1d), intent(inout) :: patch
    logical, intent(out) :: ok

    integer :: offset, count

    count = sparse_patch_value_count(patch)
    ok = size(payload) == count
    if (.not. ok) return
    offset = 0
    call consume_sparse_payload_2d(payload, offset, patch%state)
    call consume_sparse_payload_1d(payload, offset, patch%temperature)
    call consume_sparse_payload_2d( &
      payload, offset, patch%left_ghost_state)
    call consume_sparse_payload_2d( &
      payload, offset, patch%right_ghost_state)
    call consume_sparse_payload_1d( &
      payload, offset, patch%left_ghost_temperature)
    call consume_sparse_payload_1d( &
      payload, offset, patch%right_ghost_temperature)
    ok = offset == size(payload)
  end subroutine unpack_sparse_patch_values

  subroutine append_sparse_payload_2d(values, payload, offset)
    real(dp), intent(in) :: values(:, :)
    real(dp), intent(inout) :: payload(:)
    integer, intent(inout) :: offset

    integer :: count

    count = size(values)
    payload(offset + 1:offset + count) = reshape(values, [count])
    offset = offset + count
  end subroutine append_sparse_payload_2d

  subroutine append_sparse_payload_1d(values, payload, offset)
    real(dp), intent(in) :: values(:)
    real(dp), intent(inout) :: payload(:)
    integer, intent(inout) :: offset

    integer :: count

    count = size(values)
    payload(offset + 1:offset + count) = values
    offset = offset + count
  end subroutine append_sparse_payload_1d

  subroutine consume_sparse_payload_2d(payload, offset, values)
    real(dp), intent(in) :: payload(:)
    integer, intent(inout) :: offset
    real(dp), intent(out) :: values(:, :)

    integer :: count

    count = size(values)
    values = reshape(payload(offset + 1:offset + count), shape(values))
    offset = offset + count
  end subroutine consume_sparse_payload_2d

  subroutine consume_sparse_payload_1d(payload, offset, values)
    real(dp), intent(in) :: payload(:)
    integer, intent(inout) :: offset
    real(dp), intent(out) :: values(:)

    integer :: count

    count = size(values)
    values = payload(offset + 1:offset + count)
    offset = offset + count
  end subroutine consume_sparse_payload_1d

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
