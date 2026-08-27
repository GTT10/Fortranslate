module amr_eb_patch_tree_reactive_2d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use state_indices_mod, only: irho, iet
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use transport_database_mod, only: gas_transport_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_species_component, &
    reactive_mass_fraction_component, &
    reactive_conserved_to_primitive
  use reactive_2d_mod, only: advance_reactive_chemistry_2d
  use reactive_boundary_2d_mod, only: reactive_boundary_set_2d
  use reactive_eb_cfl_2d_mod, only: compute_reactive_eb_cfl_timestep_2d
  use eb_geometry_2d_mod, only: eb_geometry_2d, eb_covered_cell
  use eb_reactive_reconstruction_2d_mod, only: &
    reactive_eb_exterior_state_2d
  use eb_reactive_redistribution_2d_mod, only: &
    advance_reactive_eb_state_redistributed_2d
  use eb_reactive_transport_2d_mod, only: &
    reactive_eb_transport_timestep_2d, reactive_eb_transport_fluxes_rhs_2d
  use amr_eb_hierarchy_2d_mod, only: &
    amr_eb_patch_2d, average_down_reactive_eb_state_patch_2d
  use amr_eb_flux_register_2d_mod, only: &
    amr_eb_flux_register_2d, initialize_amr_eb_flux_register_2d, &
    accumulate_coarse_eb_fluxes_2d, accumulate_fine_eb_fluxes_2d, &
    reflux_reactive_eb_state_patch_2d
  use amr_eb_reactive_2d_mod, only: &
    prolong_reactive_eb_patch_pcm_2d, &
    build_reactive_eb_patch_exterior_2d, advance_reactive_eb_level_2d
  use amr_eb_transport_2d_mod, only: recover_transport_temperature_2d
  use amr_eb_regrid_2d_mod, only: &
    amr_eb_tagging_criteria_2d, amr_eb_regrid_plan_collection_2d, &
    plan_reactive_eb_temperature_regrid_collection_2d
  use amr_eb_patch_tree_2d_mod, only: &
    amr_eb_patch_tree_level_plan_2d, amr_eb_patch_tree_topology_2d, &
    initialize_amr_eb_patch_tree_topology_2d, &
    rebuild_amr_eb_patch_tree_topology_2d
  implicit none
  private

  real(dp), parameter :: geometry_tolerance = &
    5.0e3_dp * epsilon(1.0_dp)
  real(dp), parameter :: conservation_tolerance = &
    5.0e4_dp * epsilon(1.0_dp)
  character(len=*), parameter :: patch_tree_checkpoint_magic = &
    "PELEF_REACTIVE_AMR_EB_PATCH_TREE_2D"
  integer, parameter :: patch_tree_checkpoint_schema = 1
  integer, parameter :: checkpoint_maximum_levels = 64
  integer, parameter :: checkpoint_maximum_patches = 1000000
  integer, parameter :: checkpoint_maximum_geometry_cells = 100000000

  type, public :: reactive_amr_eb_patch_tree_node_2d
    real(dp), allocatable :: state(:, :, :)
    real(dp), allocatable :: temperature(:, :)
  end type reactive_amr_eb_patch_tree_node_2d

  type, public :: reactive_amr_eb_patch_tree_level_2d
    type(reactive_amr_eb_patch_tree_node_2d), allocatable :: patches(:)
  contains
    procedure :: patch_count => reactive_amr_eb_patch_tree_level_patch_count
  end type reactive_amr_eb_patch_tree_level_2d

  type, public :: reactive_amr_eb_patch_tree_2d
    integer :: nvar = 0
    type(amr_eb_patch_tree_topology_2d) :: topology
    type(reactive_amr_eb_patch_tree_level_2d), allocatable :: levels(:)
  contains
    procedure :: level_count => reactive_amr_eb_patch_tree_level_count
    procedure :: level_patch_count => &
      reactive_amr_eb_patch_tree_level_patch_count_at
    procedure :: is_valid => reactive_amr_eb_patch_tree_is_valid
  end type reactive_amr_eb_patch_tree_2d

  abstract interface
    subroutine reactive_amr_eb_tree_geometry_builder_2d( &
        parent_geometry, coarse_i_lower, coarse_i_upper, coarse_j_lower, &
        coarse_j_upper, refinement_ratio, child_geometry, ok)
      import :: eb_geometry_2d
      type(eb_geometry_2d), intent(in) :: parent_geometry
      integer, intent(in) :: coarse_i_lower, coarse_i_upper
      integer, intent(in) :: coarse_j_lower, coarse_j_upper
      integer, intent(in) :: refinement_ratio
      type(eb_geometry_2d), intent(out) :: child_geometry
      logical, intent(out) :: ok
    end subroutine reactive_amr_eb_tree_geometry_builder_2d
  end interface

  public :: initialize_reactive_amr_eb_patch_tree_2d
  public :: synchronize_reactive_amr_eb_patch_tree_2d
  public :: rebuild_reactive_amr_eb_patch_tree_2d
  public :: plan_tagged_reactive_amr_eb_patch_tree_2d
  public :: regrid_tagged_reactive_amr_eb_patch_tree_2d
  public :: write_reactive_amr_eb_patch_tree_2d_checkpoint
  public :: read_reactive_amr_eb_patch_tree_2d_checkpoint
  public :: write_reactive_amr_eb_patch_tree_2d_csv
  public :: compute_reactive_amr_eb_patch_tree_cfl_timestep_2d
  public :: compute_reactive_amr_eb_patch_tree_timestep_2d
  public :: advance_reactive_amr_eb_patch_tree_chemistry_2d
  public :: advance_reactive_amr_eb_patch_tree_hydro_2d
  public :: advance_reactive_amr_eb_patch_tree_strang_2d
  public :: advance_reactive_amr_eb_patch_tree_full_physics_2d
  public :: advance_reactive_amr_eb_patch_tree_to_time_2d
  public :: advance_reactive_amr_eb_patch_tree_transport_euler_2d
  public :: advance_reactive_amr_eb_patch_tree_transport_2d
  public :: composite_integral_reactive_amr_eb_patch_tree_2d
  public :: composite_reactive_amr_eb_patch_subtree_integral_2d

contains

  pure integer function reactive_amr_eb_patch_tree_level_patch_count(self) &
      result(count)
    class(reactive_amr_eb_patch_tree_level_2d), intent(in) :: self

    count = 0
    if (allocated(self%patches)) count = size(self%patches)
  end function reactive_amr_eb_patch_tree_level_patch_count

  pure integer function reactive_amr_eb_patch_tree_level_count(self) &
      result(count)
    class(reactive_amr_eb_patch_tree_2d), intent(in) :: self

    count = 0
    if (allocated(self%levels)) count = size(self%levels)
  end function reactive_amr_eb_patch_tree_level_count

  pure integer function reactive_amr_eb_patch_tree_level_patch_count_at( &
      self, level) result(count)
    class(reactive_amr_eb_patch_tree_2d), intent(in) :: self
    integer, intent(in) :: level

    count = 0
    if (.not. allocated(self%levels)) return
    if (level < 0 .or. level >= size(self%levels)) return
    count = self%levels(level + 1)%patch_count()
  end function reactive_amr_eb_patch_tree_level_patch_count_at

  logical function reactive_amr_eb_patch_tree_is_valid(self) result(valid)
    class(reactive_amr_eb_patch_tree_2d), intent(in) :: self

    type(eb_geometry_2d) :: geometry
    integer :: level, patch

    valid = self%nvar >= 1 .and. self%topology%is_valid() .and. &
      allocated(self%levels)
    if (.not. valid) return
    valid = size(self%levels) == self%topology%level_count()
    if (.not. valid) return

    do level = 1, size(self%levels)
      valid = allocated(self%levels(level)%patches) .and. &
        self%levels(level)%patch_count() == &
          self%topology%level_patch_count(level - 1)
      if (.not. valid) return
      do patch = 1, self%levels(level)%patch_count()
        call patch_geometry_at(self%topology, level, patch, geometry, valid)
        if (.not. valid) return
        valid = allocated(self%levels(level)%patches(patch)%state) .and. &
          allocated(self%levels(level)%patches(patch)%temperature)
        if (.not. valid) return
        valid = all(shape(self%levels(level)%patches(patch)%state) == &
            [self%nvar, geometry%nx, geometry%ny]) .and. &
          all(shape(self%levels(level)%patches(patch)%temperature) == &
            [geometry%nx, geometry%ny])
        if (.not. valid) return
        valid = all(ieee_is_finite( &
            self%levels(level)%patches(patch)%state)) .and. &
          all(ieee_is_finite( &
            self%levels(level)%patches(patch)%temperature)) .and. &
          all(self%levels(level)%patches(patch)%temperature > 0.0_dp)
        if (.not. valid) return
      end do
    end do
  end function reactive_amr_eb_patch_tree_is_valid

  subroutine initialize_reactive_amr_eb_patch_tree_2d( &
      species, root_state, root_temperature, topology, solution, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: root_state(:, :, :), root_temperature(:, :)
    type(amr_eb_patch_tree_topology_2d), intent(in) :: topology
    type(reactive_amr_eb_patch_tree_2d), intent(out) :: solution
    logical, intent(out) :: ok

    type(reactive_amr_eb_patch_tree_2d) :: candidate
    type(eb_geometry_2d) :: parent_geometry
    logical :: local_ok
    integer :: child, nvar, parent, relation

    solution = reactive_amr_eb_patch_tree_2d()
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar < 1 .or. .not. topology%is_valid()) return
    if (any(shape(root_state) /= &
          [nvar, topology%root_geometry%nx, topology%root_geometry%ny]) &
        .or. any(shape(root_temperature) /= &
          [topology%root_geometry%nx, topology%root_geometry%ny]) .or. &
        any(.not. ieee_is_finite(root_state)) .or. &
        any(.not. ieee_is_finite(root_temperature)) .or. &
        any(root_temperature <= 0.0_dp)) return

    candidate%nvar = nvar
    candidate%topology = topology
    allocate(candidate%levels(topology%level_count()))
    allocate(candidate%levels(1)%patches(1))
    allocate(candidate%levels(1)%patches(1)%state, source=root_state)
    allocate(candidate%levels(1)%patches(1)%temperature, &
      source=root_temperature)
    call recover_patch_temperature( &
      species, candidate%levels(1)%patches(1), topology%root_geometry, &
      local_ok)
    if (.not. local_ok) return

    do relation = 1, size(topology%relations)
      allocate(candidate%levels(relation + 1)%patches( &
        topology%relations(relation)%child_patch_count()))
      do child = 1, topology%relations(relation)%child_patch_count()
        parent = topology%relations(relation)%children(child)%parent_patch
        call patch_geometry_at( &
          topology, relation, parent, parent_geometry, local_ok)
        if (.not. local_ok) return
        allocate(candidate%levels(relation + 1)%patches(child)%state( &
          nvar, topology%relations(relation)%children(child)%geometry%nx, &
          topology%relations(relation)%children(child)%geometry%ny))
        allocate(candidate%levels(relation + 1)%patches(child)%temperature( &
          topology%relations(relation)%children(child)%geometry%nx, &
          topology%relations(relation)%children(child)%geometry%ny))
        call prolong_reactive_eb_patch_pcm_2d( &
          species, candidate%levels(relation)%patches(parent)%state, &
          candidate%levels(relation)%patches(parent)%temperature, &
          parent_geometry, &
          topology%relations(relation)%children(child)%geometry, &
          topology%relations(relation)%children(child)%patch, &
          candidate%levels(relation + 1)%patches(child)%state, &
          candidate%levels(relation + 1)%patches(child)%temperature, &
          local_ok)
        if (.not. local_ok) return
      end do
    end do

    if (.not. candidate%is_valid()) return
    solution = candidate
    ok = .true.
  end subroutine initialize_reactive_amr_eb_patch_tree_2d

  subroutine synchronize_reactive_amr_eb_patch_tree_2d( &
      species, solution, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_amr_eb_patch_tree_2d), intent(inout) :: solution
    logical, intent(out) :: ok

    type(reactive_amr_eb_patch_tree_2d) :: candidate

    ok = .false.
    if (.not. solution%is_valid()) return
    if (solution%nvar /= reactive_nvar(size(species))) return
    candidate = solution
    call synchronize_candidate(species, candidate, ok)
    if (.not. ok) return
    if (.not. candidate%is_valid()) then
      ok = .false.
      return
    end if
    solution = candidate
  end subroutine synchronize_reactive_amr_eb_patch_tree_2d

  subroutine plan_tagged_reactive_amr_eb_patch_tree_2d( &
      species, solution, criteria, maximum_levels, refinement_ratio, &
      geometry_builder, plans, tagged_cells, ok, failure_context)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_amr_eb_patch_tree_2d), intent(in) :: solution
    type(amr_eb_tagging_criteria_2d), intent(in) :: criteria
    integer, intent(in) :: maximum_levels, refinement_ratio
    procedure(reactive_amr_eb_tree_geometry_builder_2d) :: geometry_builder
    type(amr_eb_patch_tree_level_plan_2d), allocatable, intent(out) :: plans(:)
    integer, intent(out) :: tagged_cells
    logical, intent(out) :: ok
    character(len=*), intent(out), optional :: failure_context

    type(reactive_amr_eb_patch_tree_2d) :: candidate, next_candidate
    type(amr_eb_patch_tree_topology_2d) :: candidate_topology
    type(amr_eb_patch_tree_level_plan_2d), allocatable :: workspace(:)
    type(amr_eb_regrid_plan_collection_2d), allocatable :: collections(:)
    type(eb_geometry_2d) :: parent_geometry
    logical, allocatable :: tags(:, :)
    integer :: child, child_count, entry, parent, parent_count
    integer :: relation, relation_count
    logical :: local_ok

    ok = .false.
    tagged_cells = 0
    if (present(failure_context)) failure_context = "input validation"
    if (.not. solution%is_valid() .or. &
        solution%nvar /= reactive_nvar(size(species)) .or. &
        maximum_levels < 1 .or. refinement_ratio < 2 .or. &
        .not. criteria%is_valid( &
          solution%topology%root_geometry%nx, &
          solution%topology%root_geometry%ny)) return

    candidate = solution
    if (present(failure_context)) failure_context = "source synchronization"
    call synchronize_candidate(species, candidate, local_ok)
    if (.not. local_ok .or. .not. candidate%is_valid()) return

    allocate(workspace(maximum_levels - 1))
    relation_count = 0
    do relation = 1, maximum_levels - 1
      parent_count = candidate%levels(relation)%patch_count()
      if (parent_count < 1) return
      allocate(collections(parent_count))
      child_count = 0
      do parent = 1, parent_count
        if (present(failure_context)) &
          write(failure_context, '(a,i0,a,i0)') &
            "tagging level ", relation - 1, " patch ", parent
        call patch_geometry_at( &
          candidate%topology, relation, parent, parent_geometry, local_ok)
        if (.not. local_ok) return
        if (.not. criteria%is_valid( &
            parent_geometry%nx, parent_geometry%ny)) cycle
        allocate(tags(parent_geometry%nx, parent_geometry%ny))
        call plan_reactive_eb_temperature_regrid_collection_2d( &
          candidate%levels(relation)%patches(parent)%temperature, &
          parent_geometry, criteria, tags, collections(parent), local_ok)
        deallocate(tags)
        if (.not. local_ok) return
        tagged_cells = tagged_cells + collections(parent)%tagged_cell_count
        child_count = child_count + collections(parent)%patch_count()
      end do
      if (child_count == 0) then
        deallocate(collections)
        exit
      end if

      workspace(relation)%refinement_ratio = refinement_ratio
      allocate(workspace(relation)%children(child_count))
      entry = 0
      do parent = 1, parent_count
        call patch_geometry_at( &
          candidate%topology, relation, parent, parent_geometry, local_ok)
        if (.not. local_ok) return
        do child = 1, collections(parent)%patch_count()
          entry = entry + 1
          workspace(relation)%children(entry)%parent_patch = parent
          workspace(relation)%children(entry)%coarse_i_lower = &
            collections(parent)%plans(child)%coarse_i_lower
          workspace(relation)%children(entry)%coarse_i_upper = &
            collections(parent)%plans(child)%coarse_i_upper
          workspace(relation)%children(entry)%coarse_j_lower = &
            collections(parent)%plans(child)%coarse_j_lower
          workspace(relation)%children(entry)%coarse_j_upper = &
            collections(parent)%plans(child)%coarse_j_upper
          if (present(failure_context)) &
            write(failure_context, '(a,i0,a,i0)') &
              "geometry level ", relation, " patch ", entry
          call geometry_builder( &
            parent_geometry, &
            workspace(relation)%children(entry)%coarse_i_lower, &
            workspace(relation)%children(entry)%coarse_i_upper, &
            workspace(relation)%children(entry)%coarse_j_lower, &
            workspace(relation)%children(entry)%coarse_j_upper, &
            refinement_ratio, &
            workspace(relation)%children(entry)%geometry, local_ok)
          if (.not. local_ok) return
        end do
      end do
      deallocate(collections)
      relation_count = relation

      if (present(failure_context)) failure_context = "candidate topology"
      call initialize_amr_eb_patch_tree_topology_2d( &
        solution%topology%root_geometry, workspace(1:relation_count), &
        candidate_topology, local_ok)
      if (.not. local_ok) return
      if (present(failure_context)) failure_context = "candidate fields"
      call initialize_reactive_amr_eb_patch_tree_2d( &
        species, candidate%levels(1)%patches(1)%state, &
        candidate%levels(1)%patches(1)%temperature, candidate_topology, &
        next_candidate, local_ok)
      if (.not. local_ok) return
      candidate = next_candidate
    end do

    allocate(plans(relation_count))
    if (relation_count > 0) plans = workspace(1:relation_count)
    ok = .true.
    if (present(failure_context)) failure_context = "none"
  end subroutine plan_tagged_reactive_amr_eb_patch_tree_2d

  subroutine regrid_tagged_reactive_amr_eb_patch_tree_2d( &
      species, solution, criteria, maximum_levels, refinement_ratio, &
      geometry_builder, ok, changed, tagged_cells, failure_context)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_amr_eb_patch_tree_2d), intent(inout) :: solution
    type(amr_eb_tagging_criteria_2d), intent(in) :: criteria
    integer, intent(in) :: maximum_levels, refinement_ratio
    procedure(reactive_amr_eb_tree_geometry_builder_2d) :: geometry_builder
    logical, intent(out) :: ok, changed
    integer, intent(out) :: tagged_cells
    character(len=*), intent(out), optional :: failure_context

    type(amr_eb_patch_tree_level_plan_2d), allocatable :: plans(:)
    character(len=160) :: context

    ok = .false.
    changed = .false.
    tagged_cells = 0
    context = "tag planning"
    if (present(failure_context)) failure_context = context
    call plan_tagged_reactive_amr_eb_patch_tree_2d( &
      species, solution, criteria, maximum_levels, refinement_ratio, &
      geometry_builder, plans, tagged_cells, ok, context)
    if (.not. ok) then
      if (present(failure_context)) failure_context = context
      return
    end if

    context = "tree rebuild"
    call rebuild_reactive_amr_eb_patch_tree_2d( &
      species, solution, plans, ok, changed, context)
    if (.not. ok) then
      changed = .false.
      if (present(failure_context)) failure_context = context
      return
    end if
    if (present(failure_context)) failure_context = "none"
  end subroutine regrid_tagged_reactive_amr_eb_patch_tree_2d

  subroutine write_reactive_amr_eb_patch_tree_2d_checkpoint( &
      path, species, solution, time, steps, regrids, minimum_dt, ok)
    character(len=*), intent(in) :: path
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_amr_eb_patch_tree_2d), intent(in) :: solution
    real(dp), intent(in) :: time, minimum_dt
    integer, intent(in) :: steps, regrids
    logical, intent(out) :: ok

    integer :: unit, status, species_index, relation, child
    integer :: level, patch, i, j

    ok = .false.
    if (len_trim(path) == 0 .or. size(species) < 1 .or. &
        .not. solution%is_valid() .or. &
        solution%nvar /= reactive_nvar(size(species)) .or. &
        solution%level_count() > checkpoint_maximum_levels .or. &
        .not. patch_tree_checkpoint_metadata_is_valid( &
          time, steps, regrids, minimum_dt)) return

    open(newunit=unit, file=trim(path), status="replace", action="write", &
      form="formatted", iostat=status)
    if (status /= 0) return
    write(unit, '(a)', iostat=status) patch_tree_checkpoint_magic
    if (status /= 0) go to 900
    write(unit, '(*(i0,1x))', iostat=status) &
      patch_tree_checkpoint_schema, size(species), solution%nvar, &
      solution%level_count()
    if (status /= 0) go to 900
    do species_index = 1, size(species)
      write(unit, '(a)', iostat=status) trim(species(species_index)%name)
      if (status /= 0) go to 900
    end do

    call write_patch_tree_checkpoint_geometry_2d( &
      unit, solution%topology%root_geometry, status)
    if (status /= 0) go to 900
    do relation = 1, size(solution%topology%relations)
      write(unit, '(*(i0,1x))', iostat=status) relation, &
        solution%topology%relations(relation)%refinement_ratio, &
        solution%topology%relations(relation)%child_patch_count()
      if (status /= 0) go to 900
      do child = 1, &
          solution%topology%relations(relation)%child_patch_count()
        write(unit, '(*(i0,1x))', iostat=status) &
          solution%topology%relations(relation)%children(child)% &
            parent_patch, &
          solution%topology%relations(relation)%children(child)%patch% &
            coarse_i_lower, &
          solution%topology%relations(relation)%children(child)%patch% &
            coarse_i_upper, &
          solution%topology%relations(relation)%children(child)%patch% &
            coarse_j_lower, &
          solution%topology%relations(relation)%children(child)%patch% &
            coarse_j_upper
        if (status /= 0) go to 900
        call write_patch_tree_checkpoint_geometry_2d( &
          unit, solution%topology%relations(relation)%children(child)% &
            geometry, status)
        if (status /= 0) go to 900
      end do
    end do

    write(unit, '(2(es27.18e3,1x),2(i0,1x))', iostat=status) &
      time, minimum_dt, steps, regrids
    if (status /= 0) go to 900
    do level = 1, solution%level_count()
      do patch = 1, solution%levels(level)%patch_count()
        write(unit, '(*(i0,1x))', iostat=status) &
          level, patch, &
          size(solution%levels(level)%patches(patch)%state, 2), &
          size(solution%levels(level)%patches(patch)%state, 3)
        if (status /= 0) go to 900
        do j = 1, size(solution%levels(level)%patches(patch)%state, 3)
          do i = 1, size(solution%levels(level)%patches(patch)%state, 2)
            write(unit, '(*(es27.18e3,1x))', iostat=status) &
              solution%levels(level)%patches(patch)%state(:, i, j), &
              solution%levels(level)%patches(patch)%temperature(i, j)
            if (status /= 0) go to 900
          end do
        end do
      end do
    end do
    write(unit, '(a)', iostat=status) "END_CHECKPOINT"
    if (status /= 0) go to 900
    close(unit, iostat=status)
    ok = status == 0
    return

900 continue
    close(unit)
  end subroutine write_reactive_amr_eb_patch_tree_2d_checkpoint

  subroutine read_reactive_amr_eb_patch_tree_2d_checkpoint( &
      path, species, maximum_levels, solution, time, steps, regrids, &
      minimum_dt, ok)
    character(len=*), intent(in) :: path
    type(nasa7_species), intent(in) :: species(:)
    integer, intent(in) :: maximum_levels
    type(reactive_amr_eb_patch_tree_2d), intent(out) :: solution
    real(dp), intent(out) :: time, minimum_dt
    integer, intent(out) :: steps, regrids
    logical, intent(out) :: ok

    type(amr_eb_patch_tree_level_plan_2d), allocatable :: plans(:)
    type(amr_eb_patch_tree_topology_2d) :: topology
    type(reactive_amr_eb_patch_tree_2d) :: candidate
    type(eb_geometry_2d) :: root_geometry, geometry
    character(len=1024) :: magic, stored_name, end_marker
    logical :: local_ok
    integer :: unit, status, schema, stored_species, stored_nvar
    integer :: stored_levels, stored_relation, relation_patches
    integer :: species_index, relation, child, level, patch, i, j
    integer :: stored_level, stored_patch, stored_nx, stored_ny

    solution = reactive_amr_eb_patch_tree_2d()
    time = 0.0_dp
    minimum_dt = 0.0_dp
    steps = 0
    regrids = 0
    ok = .false.
    if (len_trim(path) == 0 .or. size(species) < 1 .or. &
        maximum_levels < 1 .or. &
        maximum_levels > checkpoint_maximum_levels) return

    open(newunit=unit, file=trim(path), status="old", action="read", &
      form="formatted", iostat=status)
    if (status /= 0) return
    read(unit, '(a)', iostat=status) magic
    if (status /= 0 .or. trim(magic) /= patch_tree_checkpoint_magic) &
      go to 900
    read(unit, *, iostat=status) schema, stored_species, stored_nvar, &
      stored_levels
    if (status /= 0 .or. schema /= patch_tree_checkpoint_schema .or. &
        stored_species /= size(species) .or. &
        stored_nvar /= reactive_nvar(size(species)) .or. &
        stored_levels < 1 .or. stored_levels > maximum_levels) go to 900
    do species_index = 1, stored_species
      read(unit, '(a)', iostat=status) stored_name
      if (status /= 0 .or. &
          trim(stored_name) /= trim(species(species_index)%name)) go to 900
    end do

    call read_patch_tree_checkpoint_geometry_2d( &
      unit, root_geometry, status)
    if (status /= 0) go to 900
    allocate(plans(stored_levels - 1))
    do relation = 1, size(plans)
      read(unit, *, iostat=status) stored_relation, &
        plans(relation)%refinement_ratio, relation_patches
      if (status /= 0 .or. stored_relation /= relation .or. &
          plans(relation)%refinement_ratio < 2 .or. &
          relation_patches < 1 .or. &
          relation_patches > checkpoint_maximum_patches) go to 900
      allocate(plans(relation)%children(relation_patches))
      do child = 1, relation_patches
        read(unit, *, iostat=status) &
          plans(relation)%children(child)%parent_patch, &
          plans(relation)%children(child)%coarse_i_lower, &
          plans(relation)%children(child)%coarse_i_upper, &
          plans(relation)%children(child)%coarse_j_lower, &
          plans(relation)%children(child)%coarse_j_upper
        if (status /= 0) go to 900
        call read_patch_tree_checkpoint_geometry_2d( &
          unit, plans(relation)%children(child)%geometry, status)
        if (status /= 0) go to 900
      end do
    end do
    call initialize_amr_eb_patch_tree_topology_2d( &
      root_geometry, plans, topology, local_ok)
    if (.not. local_ok) go to 900

    read(unit, *, iostat=status) time, minimum_dt, steps, regrids
    if (status /= 0 .or. .not. patch_tree_checkpoint_metadata_is_valid( &
        time, steps, regrids, minimum_dt)) go to 900

    candidate%nvar = stored_nvar
    candidate%topology = topology
    allocate(candidate%levels(stored_levels))
    do level = 1, stored_levels
      allocate(candidate%levels(level)%patches( &
        topology%level_patch_count(level - 1)))
      do patch = 1, candidate%levels(level)%patch_count()
        call patch_geometry_at( &
          topology, level, patch, geometry, local_ok)
        if (.not. local_ok) go to 900
        read(unit, *, iostat=status) stored_level, stored_patch, &
          stored_nx, stored_ny
        if (status /= 0 .or. stored_level /= level .or. &
            stored_patch /= patch .or. stored_nx /= geometry%nx .or. &
            stored_ny /= geometry%ny) go to 900
        allocate(candidate%levels(level)%patches(patch)%state( &
          stored_nvar, geometry%nx, geometry%ny))
        allocate(candidate%levels(level)%patches(patch)%temperature( &
          geometry%nx, geometry%ny))
        do j = 1, geometry%ny
          do i = 1, geometry%nx
            read(unit, *, iostat=status) &
              candidate%levels(level)%patches(patch)%state(:, i, j), &
              candidate%levels(level)%patches(patch)%temperature(i, j)
            if (status /= 0) go to 900
          end do
        end do
        if (any(.not. ieee_is_finite( &
              candidate%levels(level)%patches(patch)%state)) .or. &
            any(.not. ieee_is_finite( &
              candidate%levels(level)%patches(patch)%temperature)) .or. &
            any(candidate%levels(level)%patches(patch)%temperature <= &
              0.0_dp)) go to 900
        call recover_patch_temperature( &
          species, candidate%levels(level)%patches(patch), geometry, local_ok)
        if (.not. local_ok) go to 900
      end do
    end do
    read(unit, '(a)', iostat=status) end_marker
    if (status /= 0 .or. trim(end_marker) /= "END_CHECKPOINT") go to 900
    close(unit, iostat=status)
    if (status /= 0 .or. .not. candidate%is_valid()) then
      solution = reactive_amr_eb_patch_tree_2d()
      time = 0.0_dp
      minimum_dt = 0.0_dp
      steps = 0
      regrids = 0
      return
    end if

    solution = candidate
    ok = .true.
    return

900 continue
    close(unit)
    solution = reactive_amr_eb_patch_tree_2d()
    time = 0.0_dp
    minimum_dt = 0.0_dp
    steps = 0
    regrids = 0
  end subroutine read_reactive_amr_eb_patch_tree_2d_checkpoint

  subroutine write_reactive_amr_eb_patch_tree_2d_csv( &
      path, species, solution, time, ok)
    character(len=*), intent(in) :: path
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_amr_eb_patch_tree_2d), intent(in) :: solution
    real(dp), intent(in) :: time
    logical, intent(out) :: ok

    type(eb_geometry_2d) :: geometry
    real(dp), allocatable :: primitive(:)
    logical, allocatable :: refined(:, :)
    real(dp) :: local_temperature, sound_speed, x, y
    logical :: local_ok
    integer :: child, i, j, k, level, patch, relation, status, unit

    ok = .false.
    if (len_trim(path) == 0 .or. size(species) < 1 .or. &
        .not. solution%is_valid() .or. &
        solution%nvar /= reactive_nvar(size(species)) .or. &
        .not. ieee_is_finite(time) .or. time < 0.0_dp) return
    allocate(primitive(reactive_nprim(size(species))))
    open(newunit=unit, file=trim(path), status="replace", action="write", &
      form="formatted", iostat=status)
    if (status /= 0) return
    write(unit, '(a)', advance='no', iostat=status) &
      "level,patch,i,j,cell_dx,cell_dy,time,x,y,volume_fraction," // &
      "cell_type,boundary_length,boundary_normal_x,boundary_normal_y," // &
      "rho,u,v,w,pressure,temperature,rhoE"
    if (status /= 0) go to 900
    do k = 1, size(species)
      write(unit, '(a)', advance='no', iostat=status) &
        ",Y_" // trim(species(k)%name)
      if (status /= 0) go to 900
    end do
    write(unit, '(a)', iostat=status) ""
    if (status /= 0) go to 900

    do level = 1, solution%level_count()
      do patch = 1, solution%levels(level)%patch_count()
        call patch_geometry_at( &
          solution%topology, level, patch, geometry, local_ok)
        if (.not. local_ok) go to 900
        allocate(refined(geometry%nx, geometry%ny), source=.false.)
        if (level < solution%level_count()) then
          relation = level
          do child = 1, &
              solution%topology%relations(relation)%child_patch_count()
            if (solution%topology%relations(relation)%children(child)% &
                parent_patch /= patch) cycle
            refined( &
              solution%topology%relations(relation)%children(child)%patch% &
                coarse_i_lower: &
              solution%topology%relations(relation)%children(child)%patch% &
                coarse_i_upper, &
              solution%topology%relations(relation)%children(child)%patch% &
                coarse_j_lower: &
              solution%topology%relations(relation)%children(child)%patch% &
                coarse_j_upper) = .true.
          end do
        end if
        do j = 1, geometry%ny
          y = geometry%y_lower + (real(j, dp) - 0.5_dp) * geometry%dy
          do i = 1, geometry%nx
            if (refined(i, j)) cycle
            x = geometry%x_lower + (real(i, dp) - 0.5_dp) * geometry%dx
            call reactive_conserved_to_primitive( &
              species, solution%levels(level)%patches(patch)%state(:, i, j), &
              solution%levels(level)%patches(patch)%temperature(i, j), &
              primitive, local_temperature, sound_speed, local_ok)
            if (.not. local_ok) go to 900
            write(unit, '(*(g0,:,","))', iostat=status) &
              level - 1, patch, i, j, geometry%dx, geometry%dy, time, &
              x, y, geometry%volume_fraction(i, j), &
              geometry%cell_type(i, j), geometry%boundary_length(i, j), &
              geometry%boundary_normal_x(i, j), &
              geometry%boundary_normal_y(i, j), &
              solution%levels(level)%patches(patch)%state(irho, i, j), &
              primitive(2), primitive(3), primitive(4), primitive(5), &
              local_temperature, &
              solution%levels(level)%patches(patch)%state(iet, i, j), &
              (primitive(reactive_mass_fraction_component(k)), &
                k = 1, size(species))
            if (status /= 0) go to 900
          end do
        end do
        deallocate(refined)
      end do
    end do
    close(unit, iostat=status)
    ok = status == 0
    return

900 continue
    close(unit)
  end subroutine write_reactive_amr_eb_patch_tree_2d_csv

  subroutine rebuild_reactive_amr_eb_patch_tree_2d( &
      species, solution, plans, ok, changed, failure_context)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_amr_eb_patch_tree_2d), intent(inout) :: solution
    type(amr_eb_patch_tree_level_plan_2d), intent(in) :: plans(:)
    logical, intent(out) :: ok, changed
    character(len=*), intent(out), optional :: failure_context

    type(reactive_amr_eb_patch_tree_2d) :: collapsed, candidate
    type(amr_eb_patch_tree_topology_2d) :: new_topology
    type(eb_geometry_2d) :: parent_geometry
    real(dp), allocatable :: old_integral(:), new_integral(:)
    real(dp) :: integral_scale
    logical :: local_ok, topology_changed
    integer :: child, level, old_patch, parent, relation
    logical, allocatable :: copied(:, :)

    ok = .false.
    changed = .false.
    if (present(failure_context)) failure_context = "source validation"
    if (.not. solution%is_valid()) return
    if (solution%nvar /= reactive_nvar(size(species))) return

    if (present(failure_context)) failure_context = "topology rebuild"
    new_topology = solution%topology
    call rebuild_amr_eb_patch_tree_topology_2d( &
      new_topology, plans, local_ok, topology_changed)
    if (.not. local_ok) return
    if (.not. topology_changed) then
      ok = .true.
      if (present(failure_context)) failure_context = "none"
      return
    end if

    if (present(failure_context)) failure_context = "source integral"
    allocate(old_integral(solution%nvar), new_integral(solution%nvar))
    call composite_integral_reactive_amr_eb_patch_tree_2d( &
      solution, old_integral, local_ok)
    if (.not. local_ok) return
    if (present(failure_context)) failure_context = "source synchronization"
    collapsed = solution
    call synchronize_candidate(species, collapsed, local_ok)
    if (.not. local_ok) return
    if (present(failure_context)) failure_context = "candidate initialization"
    call initialize_reactive_amr_eb_patch_tree_2d( &
      species, collapsed%levels(1)%patches(1)%state, &
      collapsed%levels(1)%patches(1)%temperature, new_topology, &
      candidate, local_ok)
    if (.not. local_ok) return

    do level = 2, candidate%level_count()
      relation = level - 1
      do child = 1, candidate%levels(level)%patch_count()
        parent = candidate%topology%relations(relation)% &
          children(child)%parent_patch
        call patch_geometry_at( &
          candidate%topology, level - 1, parent, parent_geometry, local_ok)
        if (.not. local_ok) return
        if (present(failure_context)) &
          write(failure_context, '(a,i0,a,i0)') &
            "PCM initialization level ", level - 1, " patch ", child
        call prolong_reactive_eb_patch_pcm_2d( &
          species, candidate%levels(level - 1)%patches(parent)%state, &
          candidate%levels(level - 1)%patches(parent)%temperature, &
          parent_geometry, &
          candidate%topology%relations(relation)%children(child)%geometry, &
          candidate%topology%relations(relation)%children(child)%patch, &
          candidate%levels(level)%patches(child)%state, &
          candidate%levels(level)%patches(child)%temperature, local_ok)
        if (.not. local_ok) return
        if (present(failure_context)) &
          write(failure_context, '(a,i0,a,i0)') &
            "overlap retention level ", level - 1, " patch ", child
        allocate(copied( &
          size(candidate%levels(level)%patches(child)%temperature, 1), &
          size(candidate%levels(level)%patches(child)%temperature, 2)))
        copied = .false.
        if (level <= collapsed%level_count()) then
          do old_patch = 1, collapsed%levels(level)%patch_count()
            call retain_same_resolution_overlap( &
              candidate%levels(level)%patches(child), &
              candidate%topology%relations(relation)%children(child)%geometry, &
              collapsed%levels(level)%patches(old_patch), &
              collapsed%topology%relations(relation)%children(old_patch)% &
                geometry, copied, local_ok)
            if (.not. local_ok) return
          end do
        end if
        deallocate(copied)
        if (present(failure_context)) &
          write(failure_context, '(a,i0,a,i0)') &
            "temperature recovery level ", level - 1, " patch ", child
        call recover_patch_temperature( &
          species, candidate%levels(level)%patches(child), &
          candidate%topology%relations(relation)%children(child)%geometry, &
          local_ok)
        if (.not. local_ok) return
      end do
    end do

    if (present(failure_context)) failure_context = "candidate synchronization"
    call synchronize_candidate(species, candidate, local_ok)
    if (.not. local_ok) return
    if (.not. candidate%is_valid()) return
    if (present(failure_context)) failure_context = "candidate integral"
    call composite_integral_reactive_amr_eb_patch_tree_2d( &
      candidate, new_integral, local_ok)
    if (.not. local_ok) return
    integral_scale = max(1.0_dp, maxval(abs(old_integral)))
    if (present(failure_context)) failure_context = "conservation check"
    if (maxval(abs(new_integral - old_integral)) > &
        conservation_tolerance * integral_scale) return

    solution = candidate
    ok = .true.
    changed = .true.
    if (present(failure_context)) failure_context = "none"
  end subroutine rebuild_reactive_amr_eb_patch_tree_2d

  subroutine compute_reactive_amr_eb_patch_tree_cfl_timestep_2d( &
      species, solution, cfl, dt, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_amr_eb_patch_tree_2d), intent(in) :: solution
    real(dp), intent(in) :: cfl
    real(dp), intent(out) :: dt
    logical, intent(out) :: ok

    type(eb_geometry_2d) :: geometry
    real(dp) :: level_scale, node_dt, scaled_dt
    integer :: level, patch, refinement_ratio, active_nodes
    logical :: local_ok

    dt = 0.0_dp
    ok = .false.
    if (.not. solution%is_valid() .or. &
        solution%nvar /= reactive_nvar(size(species)) .or. &
        .not. ieee_is_finite(cfl) .or. cfl <= 0.0_dp .or. &
        cfl > 1.0_dp) return

    dt = huge(1.0_dp)
    level_scale = 1.0_dp
    active_nodes = 0
    do level = 1, solution%level_count()
      if (level > 1) then
        refinement_ratio = &
          solution%topology%relations(level - 1)%refinement_ratio
        if (level_scale > huge(1.0_dp) / &
            real(refinement_ratio, dp)) then
          dt = 0.0_dp
          return
        end if
        level_scale = level_scale * real(refinement_ratio, dp)
      end if
      do patch = 1, solution%levels(level)%patch_count()
        call patch_geometry_at( &
          solution%topology, level, patch, geometry, local_ok)
        if (.not. local_ok) then
          dt = 0.0_dp
          return
        end if
        if (count(geometry%cell_type /= eb_covered_cell) == 0) cycle
        active_nodes = active_nodes + 1
        call compute_reactive_eb_cfl_timestep_2d( &
          species, solution%levels(level)%patches(patch)%state, &
          solution%levels(level)%patches(patch)%temperature, geometry, &
          cfl, node_dt, local_ok)
        if (.not. local_ok) then
          dt = 0.0_dp
          return
        end if
        if (node_dt > huge(1.0_dp) / level_scale) then
          scaled_dt = huge(1.0_dp)
        else
          scaled_dt = level_scale * node_dt
        end if
        dt = min(dt, scaled_dt)
      end do
    end do
    if (active_nodes == 0) then
      dt = 0.0_dp
      return
    end if
    ok = ieee_is_finite(dt) .and. dt > 0.0_dp
    if (.not. ok) dt = 0.0_dp
  end subroutine compute_reactive_amr_eb_patch_tree_cfl_timestep_2d

  subroutine compute_reactive_amr_eb_patch_tree_timestep_2d( &
      species, transport, solution, hydro_cfl, transport_cfl, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, dt, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(reactive_amr_eb_patch_tree_2d), intent(in) :: solution
    real(dp), intent(in) :: hydro_cfl, transport_cfl
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled
    real(dp), intent(out) :: dt
    logical, intent(out) :: ok

    type(eb_geometry_2d) :: geometry
    real(dp) :: level_scale, maximum_diffusivity, node_dt, scaled_dt
    integer :: active_nodes, level, patch, refinement_ratio
    logical :: local_ok, transport_active

    dt = 0.0_dp
    ok = .false.
    transport_active = viscosity_enabled .or. thermal_conduction_enabled .or. &
      species_diffusion_enabled
    if (.not. solution%is_valid() .or. &
        solution%nvar /= reactive_nvar(size(species)) .or. &
        .not. ieee_is_finite(hydro_cfl) .or. hydro_cfl <= 0.0_dp .or. &
        hydro_cfl > 1.0_dp .or. &
        .not. ieee_is_finite(transport_cfl) .or. &
        transport_cfl <= 0.0_dp .or. transport_cfl > 0.5_dp .or. &
        (transport_active .and. size(transport) /= size(species))) return

    dt = huge(1.0_dp)
    level_scale = 1.0_dp
    active_nodes = 0
    do level = 1, solution%level_count()
      if (level > 1) then
        refinement_ratio = &
          solution%topology%relations(level - 1)%refinement_ratio
        if (level_scale > huge(1.0_dp) / &
            real(refinement_ratio, dp)) then
          dt = 0.0_dp
          return
        end if
        level_scale = level_scale * real(refinement_ratio, dp)
      end if
      do patch = 1, solution%levels(level)%patch_count()
        call patch_geometry_at( &
          solution%topology, level, patch, geometry, local_ok)
        if (.not. local_ok) then
          dt = 0.0_dp
          return
        end if
        if (count(geometry%cell_type /= eb_covered_cell) == 0) cycle
        active_nodes = active_nodes + 1
        call compute_reactive_eb_cfl_timestep_2d( &
          species, solution%levels(level)%patches(patch)%state, &
          solution%levels(level)%patches(patch)%temperature, geometry, &
          hydro_cfl, node_dt, local_ok)
        if (.not. local_ok) then
          dt = 0.0_dp
          return
        end if
        if (node_dt > huge(1.0_dp) / level_scale) then
          scaled_dt = huge(1.0_dp)
        else
          scaled_dt = level_scale * node_dt
        end if
        dt = min(dt, scaled_dt)
        if (transport_active) then
          call reactive_eb_transport_timestep_2d( &
            species, transport, &
            solution%levels(level)%patches(patch)%state, &
            solution%levels(level)%patches(patch)%temperature, geometry, &
            transport_cfl, viscosity_enabled, thermal_conduction_enabled, &
            species_diffusion_enabled, node_dt, maximum_diffusivity, local_ok)
          if (.not. local_ok) then
            dt = 0.0_dp
            return
          end if
          if (node_dt > huge(1.0_dp) / level_scale) then
            scaled_dt = huge(1.0_dp)
          else
            scaled_dt = level_scale * node_dt
          end if
          dt = min(dt, scaled_dt)
        end if
      end do
    end do
    if (active_nodes == 0) then
      dt = 0.0_dp
      return
    end if
    ok = ieee_is_finite(dt) .and. dt > 0.0_dp .and. dt < huge(1.0_dp)
    if (.not. ok) dt = 0.0_dp
  end subroutine compute_reactive_amr_eb_patch_tree_timestep_2d

  subroutine advance_reactive_amr_eb_patch_tree_chemistry_2d( &
      species, reactions, solution, interval, rtol, atol, ok, &
      failure_context, level_advances)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(reactive_amr_eb_patch_tree_2d), intent(inout) :: solution
    real(dp), intent(in) :: interval, rtol, atol
    logical, intent(out) :: ok
    character(len=*), intent(out), optional :: failure_context
    integer, intent(out), optional :: level_advances(:)

    type(reactive_amr_eb_patch_tree_2d) :: candidate
    integer, allocatable :: candidate_advances(:)
    character(len=160) :: context
    logical :: local_ok

    ok = .false.
    context = "input validation"
    if (present(failure_context)) failure_context = context
    if (present(level_advances)) level_advances = 0
    if (.not. valid_patch_tree_chemistry_inputs( &
          species, reactions, solution, interval, rtol, atol)) return
    if (present(level_advances)) then
      if (size(level_advances) /= solution%level_count()) return
    end if

    candidate = solution
    allocate(candidate_advances(candidate%level_count()), source=0)
    call advance_reactive_amr_eb_patch_tree_chemistry_candidate_2d( &
      species, reactions, candidate, interval, rtol, atol, &
      candidate_advances, context, local_ok)
    if (.not. local_ok) then
      if (present(failure_context)) failure_context = context
      return
    end if

    context = "final hierarchy synchronization"
    call synchronize_candidate(species, candidate, local_ok)
    if (.not. local_ok .or. .not. candidate%is_valid()) then
      if (present(failure_context)) failure_context = context
      return
    end if
    solution = candidate
    if (present(level_advances)) level_advances = candidate_advances
    ok = .true.
    if (present(failure_context)) failure_context = "none"
  end subroutine advance_reactive_amr_eb_patch_tree_chemistry_2d

  subroutine advance_reactive_amr_eb_patch_tree_strang_2d( &
      species, reactions, solution, solver, reconstruction, limiter, &
      state_redist_max_order, dt, chemistry_enabled, rtol, atol, ok, &
      state_redist_target_volume_fraction, failure_context, &
      chemistry_level_advances, hydro_level_advances)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(reactive_amr_eb_patch_tree_2d), intent(inout) :: solution
    character(len=*), intent(in) :: solver, reconstruction, limiter
    integer, intent(in) :: state_redist_max_order
    real(dp), intent(in) :: dt, rtol, atol
    logical, intent(in) :: chemistry_enabled
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: state_redist_target_volume_fraction
    character(len=*), intent(out), optional :: failure_context
    integer, intent(out), optional :: chemistry_level_advances(:)
    integer, intent(out), optional :: hydro_level_advances(:)

    type(reactive_amr_eb_patch_tree_2d) :: candidate
    integer, allocatable :: candidate_chemistry_advances(:)
    integer, allocatable :: candidate_hydro_advances(:)
    character(len=160) :: context
    real(dp) :: selected_target
    logical :: local_ok

    ok = .false.
    context = "input validation"
    if (present(failure_context)) failure_context = context
    if (present(chemistry_level_advances)) chemistry_level_advances = 0
    if (present(hydro_level_advances)) hydro_level_advances = 0
    selected_target = 0.5_dp
    if (present(state_redist_target_volume_fraction)) &
      selected_target = state_redist_target_volume_fraction
    if (.not. solution%is_valid() .or. &
        solution%nvar /= reactive_nvar(size(species)) .or. &
        .not. ieee_is_finite(dt) .or. dt <= 0.0_dp .or. &
        .not. ieee_is_finite(selected_target) .or. &
        selected_target <= 0.0_dp .or. selected_target > 1.0_dp) return
    if (chemistry_enabled) then
      if (.not. valid_patch_tree_chemistry_inputs( &
            species, reactions, solution, 0.5_dp * dt, rtol, atol)) return
    end if
    if (present(chemistry_level_advances)) then
      if (size(chemistry_level_advances) /= solution%level_count()) return
    end if
    if (present(hydro_level_advances)) then
      if (size(hydro_level_advances) /= solution%level_count()) return
    end if

    candidate = solution
    allocate(candidate_chemistry_advances(candidate%level_count()), source=0)
    allocate(candidate_hydro_advances(candidate%level_count()), source=0)
    if (chemistry_enabled) then
      call advance_reactive_amr_eb_patch_tree_chemistry_candidate_2d( &
        species, reactions, candidate, 0.5_dp * dt, rtol, atol, &
        candidate_chemistry_advances, context, local_ok)
      if (.not. local_ok) then
        if (present(failure_context)) failure_context = context
        return
      end if
    end if

    call advance_reactive_amr_eb_patch_tree_hydro_2d( &
      species, candidate, solver, reconstruction, limiter, &
      state_redist_max_order, dt, local_ok, selected_target, context, &
      candidate_hydro_advances)
    if (.not. local_ok) then
      if (present(failure_context)) failure_context = context
      return
    end if

    if (chemistry_enabled) then
      call advance_reactive_amr_eb_patch_tree_chemistry_candidate_2d( &
        species, reactions, candidate, 0.5_dp * dt, rtol, atol, &
        candidate_chemistry_advances, context, local_ok)
      if (.not. local_ok) then
        if (present(failure_context)) failure_context = context
        return
      end if
      context = "final hierarchy synchronization"
      call synchronize_candidate(species, candidate, local_ok)
      if (.not. local_ok) then
        if (present(failure_context)) failure_context = context
        return
      end if
    end if

    if (.not. candidate%is_valid()) then
      if (present(failure_context)) failure_context = "final validation"
      return
    end if
    solution = candidate
    if (present(chemistry_level_advances)) &
      chemistry_level_advances = candidate_chemistry_advances
    if (present(hydro_level_advances)) &
      hydro_level_advances = candidate_hydro_advances
    ok = .true.
    if (present(failure_context)) failure_context = "none"
  end subroutine advance_reactive_amr_eb_patch_tree_strang_2d

  subroutine advance_reactive_amr_eb_patch_tree_full_physics_2d( &
      species, reactions, transport, solution, solver, reconstruction, &
      limiter, state_redist_max_order, dt, chemistry_enabled, rtol, atol, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, &
      state_redist_target_volume_fraction, minimum_transport_theta, ok, &
      failure_context, chemistry_level_advances, transport_level_advances, &
      hydro_level_advances)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(reactive_amr_eb_patch_tree_2d), intent(inout) :: solution
    character(len=*), intent(in) :: solver, reconstruction, limiter
    integer, intent(in) :: state_redist_max_order
    real(dp), intent(in) :: dt, rtol, atol
    logical, intent(in) :: chemistry_enabled
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    real(dp), intent(in) :: state_redist_target_volume_fraction
    real(dp), intent(out) :: minimum_transport_theta
    logical, intent(out) :: ok
    character(len=*), intent(out), optional :: failure_context
    integer, intent(out), optional :: chemistry_level_advances(:)
    integer, intent(out), optional :: transport_level_advances(:)
    integer, intent(out), optional :: hydro_level_advances(:)

    type(reactive_amr_eb_patch_tree_2d) :: candidate
    integer, allocatable :: candidate_chemistry_advances(:)
    integer, allocatable :: candidate_transport_advances(:)
    integer, allocatable :: candidate_hydro_advances(:)
    integer, allocatable :: stage_transport_advances(:)
    character(len=160) :: context
    real(dp) :: first_transport_theta, second_transport_theta
    logical :: local_ok, transport_active

    ok = .false.
    minimum_transport_theta = 1.0_dp
    context = "input validation"
    if (present(failure_context)) failure_context = context
    if (present(chemistry_level_advances)) chemistry_level_advances = 0
    if (present(transport_level_advances)) transport_level_advances = 0
    if (present(hydro_level_advances)) hydro_level_advances = 0
    if (.not. solution%is_valid() .or. &
        solution%nvar /= reactive_nvar(size(species)) .or. &
        .not. ieee_is_finite(dt) .or. dt <= 0.0_dp .or. &
        .not. ieee_is_finite(state_redist_target_volume_fraction) .or. &
        state_redist_target_volume_fraction <= 0.0_dp .or. &
        state_redist_target_volume_fraction > 1.0_dp) return
    if (chemistry_enabled) then
      if (.not. valid_patch_tree_chemistry_inputs( &
            species, reactions, solution, 0.5_dp * dt, rtol, atol)) return
    end if
    transport_active = viscosity_enabled .or. thermal_conduction_enabled .or. &
      species_diffusion_enabled
    if (transport_active .and. size(transport) /= size(species)) return
    if (present(chemistry_level_advances)) then
      if (size(chemistry_level_advances) /= solution%level_count()) return
    end if
    if (present(transport_level_advances)) then
      if (size(transport_level_advances) /= solution%level_count()) return
    end if
    if (present(hydro_level_advances)) then
      if (size(hydro_level_advances) /= solution%level_count()) return
    end if

    candidate = solution
    allocate(candidate_chemistry_advances(candidate%level_count()), source=0)
    allocate(candidate_transport_advances(candidate%level_count()), source=0)
    allocate(candidate_hydro_advances(candidate%level_count()), source=0)
    allocate(stage_transport_advances(candidate%level_count()), source=0)
    if (chemistry_enabled) then
      call advance_reactive_amr_eb_patch_tree_chemistry_candidate_2d( &
        species, reactions, candidate, 0.5_dp * dt, rtol, atol, &
        candidate_chemistry_advances, context, local_ok)
      if (.not. local_ok) then
        if (present(failure_context)) failure_context = context
        return
      end if
    end if

    call advance_reactive_amr_eb_patch_tree_transport_2d( &
      species, transport, candidate, 0.5_dp * dt, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, boundaries, &
      state_redist_target_volume_fraction, state_redist_max_order, &
      first_transport_theta, local_ok, context, stage_transport_advances)
    if (.not. local_ok) then
      if (present(failure_context)) failure_context = context
      return
    end if
    candidate_transport_advances = &
      candidate_transport_advances + stage_transport_advances

    call advance_reactive_amr_eb_patch_tree_hydro_2d( &
      species, candidate, solver, reconstruction, limiter, &
      state_redist_max_order, dt, local_ok, &
      state_redist_target_volume_fraction, context, &
      candidate_hydro_advances)
    if (.not. local_ok) then
      if (present(failure_context)) failure_context = context
      return
    end if

    stage_transport_advances = 0
    call advance_reactive_amr_eb_patch_tree_transport_2d( &
      species, transport, candidate, 0.5_dp * dt, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, boundaries, &
      state_redist_target_volume_fraction, state_redist_max_order, &
      second_transport_theta, local_ok, context, stage_transport_advances)
    if (.not. local_ok) then
      if (present(failure_context)) failure_context = context
      return
    end if
    candidate_transport_advances = &
      candidate_transport_advances + stage_transport_advances

    if (chemistry_enabled) then
      call advance_reactive_amr_eb_patch_tree_chemistry_candidate_2d( &
        species, reactions, candidate, 0.5_dp * dt, rtol, atol, &
        candidate_chemistry_advances, context, local_ok)
      if (.not. local_ok) then
        if (present(failure_context)) failure_context = context
        return
      end if
    end if
    context = "final hierarchy synchronization"
    call synchronize_candidate(species, candidate, local_ok)
    if (.not. local_ok .or. .not. candidate%is_valid()) then
      if (present(failure_context)) failure_context = context
      return
    end if

    solution = candidate
    minimum_transport_theta = min( &
      first_transport_theta, second_transport_theta)
    if (present(chemistry_level_advances)) &
      chemistry_level_advances = candidate_chemistry_advances
    if (present(transport_level_advances)) &
      transport_level_advances = candidate_transport_advances
    if (present(hydro_level_advances)) &
      hydro_level_advances = candidate_hydro_advances
    ok = .true.
    if (present(failure_context)) failure_context = "none"
  end subroutine advance_reactive_amr_eb_patch_tree_full_physics_2d

  subroutine advance_reactive_amr_eb_patch_tree_to_time_2d( &
      species, reactions, transport, solution, solver, reconstruction, &
      limiter, state_redist_max_order, time, final_time, steps, &
      maximum_steps, hydro_cfl, transport_cfl, chemistry_enabled, rtol, &
      atol, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, &
      state_redist_target_volume_fraction, minimum_dt, &
      minimum_transport_theta, ok, failure_context, advanced_steps, &
      chemistry_level_advances, transport_level_advances, &
      hydro_level_advances)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(reactive_amr_eb_patch_tree_2d), intent(inout) :: solution
    character(len=*), intent(in) :: solver, reconstruction, limiter
    integer, intent(in) :: state_redist_max_order
    real(dp), intent(inout) :: time
    real(dp), intent(in) :: final_time
    integer, intent(inout) :: steps
    integer, intent(in) :: maximum_steps
    real(dp), intent(in) :: hydro_cfl, transport_cfl, rtol, atol
    logical, intent(in) :: chemistry_enabled
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    real(dp), intent(in) :: state_redist_target_volume_fraction
    real(dp), intent(out) :: minimum_dt, minimum_transport_theta
    logical, intent(out) :: ok
    character(len=*), intent(out), optional :: failure_context
    integer, intent(out), optional :: advanced_steps
    integer, intent(out), optional :: chemistry_level_advances(:)
    integer, intent(out), optional :: transport_level_advances(:)
    integer, intent(out), optional :: hydro_level_advances(:)

    type(reactive_amr_eb_patch_tree_2d) :: candidate
    integer, allocatable :: accumulated_chemistry(:)
    integer, allocatable :: accumulated_transport(:)
    integer, allocatable :: accumulated_hydro(:)
    integer, allocatable :: step_chemistry(:), step_transport(:)
    integer, allocatable :: step_hydro(:)
    character(len=160) :: context
    real(dp) :: dt, remaining, step_theta, time_tolerance
    logical :: local_ok, transport_active

    ok = .false.
    minimum_dt = 0.0_dp
    minimum_transport_theta = 1.0_dp
    context = "input validation"
    if (present(failure_context)) failure_context = context
    if (present(advanced_steps)) advanced_steps = 0
    if (present(chemistry_level_advances)) chemistry_level_advances = 0
    if (present(transport_level_advances)) transport_level_advances = 0
    if (present(hydro_level_advances)) hydro_level_advances = 0
    transport_active = viscosity_enabled .or. thermal_conduction_enabled .or. &
      species_diffusion_enabled
    time_tolerance = 16.0_dp * epsilon(1.0_dp) * &
      max(tiny(1.0_dp), abs(final_time))
    if (.not. solution%is_valid() .or. &
        solution%nvar /= reactive_nvar(size(species)) .or. &
        .not. ieee_is_finite(time) .or. time < 0.0_dp .or. &
        .not. ieee_is_finite(final_time) .or. &
        final_time < time - time_tolerance .or. &
        steps < 0 .or. maximum_steps < steps .or. &
        .not. ieee_is_finite(hydro_cfl) .or. hydro_cfl <= 0.0_dp .or. &
        hydro_cfl > 1.0_dp .or. &
        .not. ieee_is_finite(transport_cfl) .or. &
        transport_cfl <= 0.0_dp .or. transport_cfl > 0.5_dp .or. &
        .not. ieee_is_finite(rtol) .or. rtol <= 0.0_dp .or. &
        .not. ieee_is_finite(atol) .or. atol <= 0.0_dp .or. &
        .not. ieee_is_finite(state_redist_target_volume_fraction) .or. &
        state_redist_target_volume_fraction <= 0.0_dp .or. &
        state_redist_target_volume_fraction > 1.0_dp .or. &
        (transport_active .and. size(transport) /= size(species))) return
    if (present(chemistry_level_advances)) then
      if (size(chemistry_level_advances) /= solution%level_count()) return
    end if
    if (present(transport_level_advances)) then
      if (size(transport_level_advances) /= solution%level_count()) return
    end if
    if (present(hydro_level_advances)) then
      if (size(hydro_level_advances) /= solution%level_count()) return
    end if

    allocate(accumulated_chemistry(solution%level_count()), source=0)
    allocate(accumulated_transport(solution%level_count()), source=0)
    allocate(accumulated_hydro(solution%level_count()), source=0)
    allocate(step_chemistry(solution%level_count()), source=0)
    allocate(step_transport(solution%level_count()), source=0)
    allocate(step_hydro(solution%level_count()), source=0)
    do
      remaining = final_time - time
      if (remaining <= time_tolerance) exit
      if (steps >= maximum_steps) then
        if (present(failure_context)) failure_context = "maximum steps"
        return
      end if
      context = "timestep"
      call compute_reactive_amr_eb_patch_tree_timestep_2d( &
        species, transport, solution, hydro_cfl, transport_cfl, &
        viscosity_enabled, thermal_conduction_enabled, &
        species_diffusion_enabled, dt, local_ok)
      if (.not. local_ok) then
        if (present(failure_context)) failure_context = context
        return
      end if
      dt = min(dt, remaining)
      if (.not. ieee_is_finite(dt) .or. dt <= 0.0_dp) then
        if (present(failure_context)) failure_context = context
        return
      end if

      candidate = solution
      step_chemistry = 0
      step_transport = 0
      step_hydro = 0
      call advance_reactive_amr_eb_patch_tree_full_physics_2d( &
        species, reactions, transport, candidate, solver, reconstruction, &
        limiter, state_redist_max_order, dt, chemistry_enabled, rtol, atol, &
        viscosity_enabled, thermal_conduction_enabled, &
        species_diffusion_enabled, barodiffusion_enabled, boundaries, &
        state_redist_target_volume_fraction, step_theta, local_ok, context, &
        step_chemistry, step_transport, step_hydro)
      if (.not. local_ok) then
        if (present(failure_context)) failure_context = context
        return
      end if
      if (any(step_chemistry > huge(steps) - accumulated_chemistry) .or. &
          any(step_transport > huge(steps) - accumulated_transport) .or. &
          any(step_hydro > huge(steps) - accumulated_hydro)) then
        if (present(failure_context)) failure_context = &
          "advance count overflow"
        return
      end if

      solution = candidate
      time = time + dt
      steps = steps + 1
      accumulated_chemistry = accumulated_chemistry + step_chemistry
      accumulated_transport = accumulated_transport + step_transport
      accumulated_hydro = accumulated_hydro + step_hydro
      if (minimum_dt == 0.0_dp) then
        minimum_dt = dt
      else
        minimum_dt = min(minimum_dt, dt)
      end if
      minimum_transport_theta = min(minimum_transport_theta, step_theta)
      if (present(advanced_steps)) advanced_steps = advanced_steps + 1
      if (present(chemistry_level_advances)) &
        chemistry_level_advances = accumulated_chemistry
      if (present(transport_level_advances)) &
        transport_level_advances = accumulated_transport
      if (present(hydro_level_advances)) &
        hydro_level_advances = accumulated_hydro
    end do

    time = final_time
    ok = .true.
    if (present(failure_context)) failure_context = "none"
  end subroutine advance_reactive_amr_eb_patch_tree_to_time_2d

  subroutine advance_reactive_amr_eb_patch_tree_chemistry_candidate_2d( &
      species, reactions, solution, interval, rtol, atol, level_advances, &
      failure_context, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(reactive_amr_eb_patch_tree_2d), intent(inout) :: solution
    real(dp), intent(in) :: interval, rtol, atol
    integer, intent(inout) :: level_advances(:)
    character(len=*), intent(inout) :: failure_context
    logical, intent(out) :: ok

    type(eb_geometry_2d) :: geometry
    logical, allocatable :: active_mask(:, :)
    logical :: local_ok
    integer :: level, patch

    ok = .false.
    if (.not. valid_patch_tree_chemistry_inputs( &
          species, reactions, solution, interval, rtol, atol) .or. &
        size(level_advances) /= solution%level_count()) return
    do level = 1, solution%level_count()
      do patch = 1, solution%levels(level)%patch_count()
        call patch_geometry_at( &
          solution%topology, level, patch, geometry, local_ok)
        if (.not. local_ok) return
        if (allocated(active_mask)) deallocate(active_mask)
        allocate(active_mask(geometry%nx, geometry%ny))
        active_mask = geometry%cell_type /= eb_covered_cell
        write(failure_context, '(a,i0,a,i0)') &
          "chemistry level ", level - 1, " patch ", patch
        call advance_reactive_chemistry_2d( &
          species, reactions, solution%levels(level)%patches(patch)%state, &
          solution%levels(level)%patches(patch)%temperature, geometry%nx, &
          geometry%ny, interval, rtol, atol, local_ok, active_mask)
        if (.not. local_ok) return
        level_advances(level) = level_advances(level) + 1
      end do
    end do
    ok = .true.
  end subroutine &
    advance_reactive_amr_eb_patch_tree_chemistry_candidate_2d

  logical function valid_patch_tree_chemistry_inputs( &
      species, reactions, solution, interval, rtol, atol) result(valid)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(reactive_amr_eb_patch_tree_2d), intent(in) :: solution
    real(dp), intent(in) :: interval, rtol, atol

    valid = solution%is_valid() .and. &
      solution%nvar == reactive_nvar(size(species)) .and. &
      size(reactions) >= 1 .and. ieee_is_finite(interval) .and. &
      interval >= 0.0_dp .and. ieee_is_finite(rtol) .and. rtol > 0.0_dp .and. &
      ieee_is_finite(atol) .and. atol > 0.0_dp
  end function valid_patch_tree_chemistry_inputs

  subroutine advance_reactive_amr_eb_patch_tree_transport_euler_2d( &
      species, transport, solution, interval, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, boundaries, target_volume_fraction, max_order, &
      minimum_theta, ok, failure_context, level_advances)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(reactive_amr_eb_patch_tree_2d), intent(inout) :: solution
    real(dp), intent(in) :: interval, target_volume_fraction
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    integer, intent(in) :: max_order
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok
    character(len=*), intent(out), optional :: failure_context
    integer, intent(out), optional :: level_advances(:)

    type(reactive_amr_eb_patch_tree_2d) :: candidate
    real(dp), allocatable :: x_flux(:, :, :), y_flux(:, :, :)
    integer, allocatable :: candidate_advances(:)
    character(len=160) :: context
    real(dp) :: candidate_theta
    logical :: local_ok

    ok = .false.
    minimum_theta = 1.0_dp
    context = "input validation"
    if (present(failure_context)) failure_context = context
    if (present(level_advances)) level_advances = 0
    if (.not. valid_patch_tree_transport_inputs( &
          species, transport, solution, interval, target_volume_fraction)) &
      return
    if (present(level_advances)) then
      if (size(level_advances) /= solution%level_count()) return
    end if

    candidate = solution
    candidate_theta = 1.0_dp
    allocate(candidate_advances(candidate%level_count()), source=0)
    call advance_reactive_amr_eb_patch_tree_transport_node_2d( &
      species, transport, candidate, 1, 1, interval, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, boundaries, target_volume_fraction, max_order, &
      x_flux, y_flux, candidate_theta, candidate_advances, context, local_ok)
    if (.not. local_ok) then
      if (present(failure_context)) failure_context = context
      return
    end if

    context = "final hierarchy synchronization"
    call synchronize_candidate(species, candidate, local_ok)
    if (.not. local_ok .or. .not. candidate%is_valid()) then
      if (present(failure_context)) failure_context = context
      return
    end if
    solution = candidate
    minimum_theta = candidate_theta
    if (present(level_advances)) level_advances = candidate_advances
    ok = .true.
    if (present(failure_context)) failure_context = "none"
  end subroutine advance_reactive_amr_eb_patch_tree_transport_euler_2d

  subroutine advance_reactive_amr_eb_patch_tree_transport_2d( &
      species, transport, solution, interval, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, boundaries, target_volume_fraction, max_order, &
      minimum_theta, ok, failure_context, level_advances)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(reactive_amr_eb_patch_tree_2d), intent(inout) :: solution
    real(dp), intent(in) :: interval, target_volume_fraction
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    integer, intent(in) :: max_order
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok
    character(len=*), intent(out), optional :: failure_context
    integer, intent(out), optional :: level_advances(:)

    type(reactive_amr_eb_patch_tree_2d) :: candidate
    type(eb_geometry_2d) :: geometry
    real(dp), allocatable :: temperature_work(:, :)
    integer, allocatable :: first_advances(:), second_advances(:)
    character(len=160) :: context
    real(dp) :: first_theta, second_theta
    logical :: local_ok, transport_active
    integer :: level, patch

    ok = .false.
    minimum_theta = 1.0_dp
    context = "input validation"
    if (present(failure_context)) failure_context = context
    if (present(level_advances)) level_advances = 0
    if (.not. solution%is_valid() .or. &
        solution%nvar /= reactive_nvar(size(species)) .or. &
        .not. ieee_is_finite(interval) .or. interval < 0.0_dp .or. &
        .not. ieee_is_finite(target_volume_fraction) .or. &
        target_volume_fraction <= 0.0_dp .or. &
        target_volume_fraction > 1.0_dp) return
    if (present(level_advances)) then
      if (size(level_advances) /= solution%level_count()) return
    end if
    transport_active = viscosity_enabled .or. thermal_conduction_enabled .or. &
      species_diffusion_enabled
    if (interval <= tiny(1.0_dp) .or. .not. transport_active) then
      ok = .true.
      if (present(failure_context)) failure_context = "none"
      return
    end if
    if (size(transport) /= size(species)) return

    candidate = solution
    allocate(first_advances(candidate%level_count()), source=0)
    allocate(second_advances(candidate%level_count()), source=0)
    call advance_reactive_amr_eb_patch_tree_transport_euler_2d( &
      species, transport, candidate, interval, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, boundaries, target_volume_fraction, max_order, &
      first_theta, local_ok, context, first_advances)
    if (.not. local_ok) then
      if (present(failure_context)) failure_context = context
      return
    end if
    call advance_reactive_amr_eb_patch_tree_transport_euler_2d( &
      species, transport, candidate, interval, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, boundaries, target_volume_fraction, max_order, &
      second_theta, local_ok, context, second_advances)
    if (.not. local_ok) then
      if (present(failure_context)) failure_context = context
      return
    end if

    do level = 1, candidate%level_count()
      do patch = 1, candidate%levels(level)%patch_count()
        call patch_geometry_at( &
          candidate%topology, level, patch, geometry, local_ok)
        if (.not. local_ok) then
          if (present(failure_context)) failure_context = "stage blend geometry"
          return
        end if
        candidate%levels(level)%patches(patch)%state = 0.5_dp * ( &
          solution%levels(level)%patches(patch)%state + &
          candidate%levels(level)%patches(patch)%state)
        if (allocated(temperature_work)) deallocate(temperature_work)
        allocate(temperature_work(geometry%nx, geometry%ny))
        call recover_transport_temperature_2d( &
          species, candidate%levels(level)%patches(patch)%state, &
          0.5_dp * (solution%levels(level)%patches(patch)%temperature + &
            candidate%levels(level)%patches(patch)%temperature), geometry, &
          temperature_work, local_ok)
        if (.not. local_ok) then
          if (present(failure_context)) failure_context = &
            "stage blend temperature recovery"
          return
        end if
        candidate%levels(level)%patches(patch)%temperature = temperature_work
      end do
    end do
    context = "final hierarchy synchronization"
    call synchronize_candidate(species, candidate, local_ok)
    if (.not. local_ok .or. .not. candidate%is_valid()) then
      if (present(failure_context)) failure_context = context
      return
    end if

    solution = candidate
    minimum_theta = min(first_theta, second_theta)
    if (present(level_advances)) &
      level_advances = first_advances + second_advances
    ok = .true.
    if (present(failure_context)) failure_context = "none"
  end subroutine advance_reactive_amr_eb_patch_tree_transport_2d

  recursive subroutine advance_reactive_amr_eb_patch_tree_transport_node_2d( &
      species, transport, solution, level, patch_index, interval, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, boundaries, &
      target_volume_fraction, max_order, x_flux, y_flux, minimum_theta, &
      level_advances, failure_context, ok, exterior)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(reactive_amr_eb_patch_tree_2d), intent(inout) :: solution
    integer, intent(in) :: level, patch_index, max_order
    real(dp), intent(in) :: interval, target_volume_fraction
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    type(reactive_boundary_set_2d), intent(in) :: boundaries
    real(dp), allocatable, intent(out) :: x_flux(:, :, :), y_flux(:, :, :)
    real(dp), intent(inout) :: minimum_theta
    integer, intent(inout) :: level_advances(:)
    character(len=*), intent(inout) :: failure_context
    logical, intent(out) :: ok
    type(reactive_eb_exterior_state_2d), intent(in), optional :: exterior

    type(eb_geometry_2d) :: geometry, child_geometry
    type(amr_eb_flux_register_2d), allocatable :: registers(:)
    type(reactive_eb_exterior_state_2d) :: child_exterior
    real(dp), allocatable :: state_start(:, :, :), temperature_start(:, :)
    real(dp), allocatable :: state_end(:, :, :), temperature_end(:, :)
    real(dp), allocatable :: rhs(:, :, :)
    real(dp), allocatable :: parent_work(:, :, :), parent_work_temperature(:, :)
    real(dp), allocatable :: child_work(:, :, :), child_work_temperature(:, :)
    real(dp), allocatable :: child_x_flux(:, :, :), child_y_flux(:, :, :)
    real(dp), allocatable :: integral_before(:)
    real(dp) :: alpha, child_interval, node_theta
    integer :: child, child_count, first_child, global_child, ratio, substep
    logical :: local_ok

    ok = .false.
    if (.not. ieee_is_finite(interval) .or. interval <= 0.0_dp .or. &
        level < 1 .or. level > solution%level_count() .or. &
        patch_index < 1 .or. &
        patch_index > solution%levels(level)%patch_count() .or. &
        size(level_advances) /= solution%level_count()) return
    call patch_geometry_at( &
      solution%topology, level, patch_index, geometry, local_ok)
    if (.not. local_ok) return

    child_count = 0
    first_child = 0
    if (level < solution%level_count()) then
      first_child = solution%topology%relations(level)% &
        child_offsets(patch_index) + 1
      child_count = solution%topology%relations(level)% &
          child_offsets(patch_index + 1) - &
        solution%topology%relations(level)%child_offsets(patch_index)
    end if
    if (child_count > 0) then
      allocate(integral_before(solution%nvar))
      call composite_reactive_amr_eb_patch_subtree_integral_2d( &
        solution, level, patch_index, integral_before, local_ok)
      if (.not. local_ok) return
    end if

    allocate(state_start, source= &
      solution%levels(level)%patches(patch_index)%state)
    allocate(temperature_start, source= &
      solution%levels(level)%patches(patch_index)%temperature)
    allocate(state_end, mold=state_start)
    allocate(temperature_end, mold=temperature_start)
    allocate(rhs, mold=state_start)
    allocate(x_flux(solution%nvar, 0:geometry%nx, geometry%ny))
    allocate(y_flux(solution%nvar, geometry%nx, 0:geometry%ny))
    write(failure_context, '(a,i0,a,i0)') &
      "transport level ", level - 1, " patch ", patch_index
    if (present(exterior)) then
      call reactive_eb_transport_fluxes_rhs_2d( &
        species, transport, state_start, temperature_start, geometry, &
        interval, viscosity_enabled, thermal_conduction_enabled, &
        species_diffusion_enabled, barodiffusion_enabled, boundaries, rhs, &
        x_flux, y_flux, node_theta, local_ok, exterior)
    else
      call reactive_eb_transport_fluxes_rhs_2d( &
        species, transport, state_start, temperature_start, geometry, &
        interval, viscosity_enabled, thermal_conduction_enabled, &
        species_diffusion_enabled, barodiffusion_enabled, boundaries, rhs, &
        x_flux, y_flux, node_theta, local_ok)
    end if
    if (.not. local_ok) return
    minimum_theta = min(minimum_theta, node_theta)
    call advance_reactive_eb_state_redistributed_2d( &
      species, state_start, temperature_start, geometry, rhs, interval, &
      state_end, temperature_end, local_ok, target_volume_fraction, max_order)
    if (.not. local_ok) return
    solution%levels(level)%patches(patch_index)%state = state_end
    solution%levels(level)%patches(patch_index)%temperature = temperature_end
    level_advances(level) = level_advances(level) + 1
    if (child_count == 0) then
      ok = .true.
      return
    end if

    ratio = solution%topology%relations(level)%refinement_ratio
    child_interval = interval / real(ratio, dp)
    allocate(registers(child_count))
    do child = 1, child_count
      global_child = first_child + child - 1
      child_geometry = solution%topology%relations(level)% &
        children(global_child)%geometry
      call initialize_amr_eb_flux_register_2d( &
        geometry, child_geometry, &
        solution%topology%relations(level)%children(global_child)%patch, &
        solution%nvar, registers(child), local_ok)
      if (.not. local_ok) return
      call accumulate_coarse_eb_fluxes_2d( &
        registers(child), geometry, child_geometry, &
        solution%topology%relations(level)%children(global_child)%patch, &
        x_flux, y_flux, interval, local_ok)
      if (.not. local_ok) return
    end do

    do substep = 1, ratio
      alpha = real(substep - 1, dp) / real(ratio, dp)
      do child = 1, child_count
        global_child = first_child + child - 1
        child_geometry = solution%topology%relations(level)% &
          children(global_child)%geometry
        call build_reactive_eb_patch_exterior_2d( &
          species, state_start, temperature_start, state_end, &
          temperature_end, geometry, child_geometry, &
          solution%topology%relations(level)%children(global_child)%patch, &
          alpha, child_exterior, local_ok, &
          solution%levels(level + 1)%patches(global_child)%state, &
          solution%levels(level + 1)%patches(global_child)%temperature)
        if (.not. local_ok) return
        call advance_reactive_amr_eb_patch_tree_transport_node_2d( &
          species, transport, solution, level + 1, global_child, &
          child_interval, viscosity_enabled, thermal_conduction_enabled, &
          species_diffusion_enabled, barodiffusion_enabled, boundaries, &
          target_volume_fraction, max_order, child_x_flux, child_y_flux, &
          minimum_theta, level_advances, failure_context, local_ok, &
          child_exterior)
        if (.not. local_ok) return
        call accumulate_fine_eb_fluxes_2d( &
          registers(child), geometry, child_geometry, &
          solution%topology%relations(level)%children(global_child)%patch, &
          child_x_flux, child_y_flux, child_interval, local_ok)
        if (.not. local_ok) return
      end do
    end do

    allocate(parent_work, mold=state_end)
    allocate(parent_work_temperature, mold=temperature_end)
    do child = 1, child_count
      global_child = first_child + child - 1
      child_geometry = solution%topology%relations(level)% &
        children(global_child)%geometry
      if (allocated(child_work)) deallocate(child_work)
      if (allocated(child_work_temperature)) deallocate(child_work_temperature)
      allocate(child_work, mold= &
        solution%levels(level + 1)%patches(global_child)%state)
      allocate(child_work_temperature, mold= &
        solution%levels(level + 1)%patches(global_child)%temperature)
      call reflux_reactive_eb_state_patch_2d( &
        species, solution%levels(level)%patches(patch_index)%state, &
        solution%levels(level)%patches(patch_index)%temperature, geometry, &
        solution%levels(level + 1)%patches(global_child)%state, &
        solution%levels(level + 1)%patches(global_child)%temperature, &
        child_geometry, &
        solution%topology%relations(level)%children(global_child)%patch, &
        registers(child), parent_work, parent_work_temperature, child_work, &
        child_work_temperature, local_ok)
      if (.not. local_ok) return
      solution%levels(level)%patches(patch_index)%state = parent_work
      solution%levels(level)%patches(patch_index)%temperature = &
        parent_work_temperature
      solution%levels(level + 1)%patches(global_child)%state = child_work
      solution%levels(level + 1)%patches(global_child)%temperature = &
        child_work_temperature
    end do

    do child = 1, child_count
      global_child = first_child + child - 1
      child_geometry = solution%topology%relations(level)% &
        children(global_child)%geometry
      call average_down_reactive_eb_state_patch_2d( &
        species, solution%levels(level)%patches(patch_index)%state, &
        solution%levels(level)%patches(patch_index)%temperature, geometry, &
        solution%levels(level + 1)%patches(global_child)%state, &
        child_geometry, &
        solution%topology%relations(level)%children(global_child)%patch, &
        parent_work, parent_work_temperature, local_ok)
      if (.not. local_ok) return
      solution%levels(level)%patches(patch_index)%state = parent_work
      solution%levels(level)%patches(patch_index)%temperature = &
        parent_work_temperature
    end do

    call close_reactive_amr_eb_patch_subtree_conservation_2d( &
      species, solution, level, patch_index, integral_before, x_flux, &
      y_flux, interval, local_ok)
    if (.not. local_ok) return
    ok = .true.
  end subroutine advance_reactive_amr_eb_patch_tree_transport_node_2d

  logical function valid_patch_tree_transport_inputs( &
      species, transport, solution, interval, target_volume_fraction) &
      result(valid)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(reactive_amr_eb_patch_tree_2d), intent(in) :: solution
    real(dp), intent(in) :: interval, target_volume_fraction

    valid = solution%is_valid() .and. &
      solution%nvar == reactive_nvar(size(species)) .and. &
      size(transport) == size(species) .and. ieee_is_finite(interval) .and. &
      interval > 0.0_dp .and. ieee_is_finite(target_volume_fraction) .and. &
      target_volume_fraction > 0.0_dp .and. &
      target_volume_fraction <= 1.0_dp
  end function valid_patch_tree_transport_inputs

  subroutine advance_reactive_amr_eb_patch_tree_hydro_2d( &
      species, solution, solver, reconstruction, limiter, &
      state_redist_max_order, dt, ok, &
      state_redist_target_volume_fraction, failure_context, level_advances)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_amr_eb_patch_tree_2d), intent(inout) :: solution
    character(len=*), intent(in) :: solver, reconstruction, limiter
    integer, intent(in) :: state_redist_max_order
    real(dp), intent(in) :: dt
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: state_redist_target_volume_fraction
    character(len=*), intent(out), optional :: failure_context
    integer, intent(out), optional :: level_advances(:)

    type(reactive_amr_eb_patch_tree_2d) :: candidate
    real(dp), allocatable :: x_flux(:, :, :), y_flux(:, :, :)
    real(dp) :: selected_target
    integer, allocatable :: candidate_advances(:)
    character(len=160) :: context
    logical :: local_ok

    ok = .false.
    context = "input validation"
    if (present(failure_context)) failure_context = context
    if (present(level_advances)) level_advances = 0
    selected_target = 0.5_dp
    if (present(state_redist_target_volume_fraction)) &
      selected_target = state_redist_target_volume_fraction
    if (.not. solution%is_valid() .or. &
        solution%nvar /= reactive_nvar(size(species)) .or. &
        .not. ieee_is_finite(dt) .or. dt <= 0.0_dp .or. &
        .not. ieee_is_finite(selected_target) .or. &
        selected_target <= 0.0_dp .or. selected_target > 1.0_dp) return
    if (present(level_advances)) then
      if (size(level_advances) /= solution%level_count()) return
    end if

    candidate = solution
    allocate(candidate_advances(candidate%level_count()), source=0)
    call advance_reactive_amr_eb_patch_tree_hydro_node_2d( &
      species, candidate, 1, 1, solver, reconstruction, limiter, &
      state_redist_max_order, selected_target, dt, x_flux, y_flux, &
      candidate_advances, context, local_ok)
    if (.not. local_ok) then
      if (present(failure_context)) failure_context = context
      return
    end if

    context = "final hierarchy synchronization"
    call synchronize_candidate(species, candidate, local_ok)
    if (.not. local_ok .or. .not. candidate%is_valid()) then
      if (present(failure_context)) failure_context = context
      return
    end if
    solution = candidate
    if (present(level_advances)) level_advances = candidate_advances
    ok = .true.
    if (present(failure_context)) failure_context = "none"
  end subroutine advance_reactive_amr_eb_patch_tree_hydro_2d

  recursive subroutine advance_reactive_amr_eb_patch_tree_hydro_node_2d( &
      species, solution, level, patch_index, solver, reconstruction, limiter, &
      state_redist_max_order, selected_target, interval, x_flux, y_flux, &
      level_advances, failure_context, ok, exterior)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_amr_eb_patch_tree_2d), intent(inout) :: solution
    integer, intent(in) :: level, patch_index, state_redist_max_order
    character(len=*), intent(in) :: solver, reconstruction, limiter
    real(dp), intent(in) :: selected_target, interval
    real(dp), allocatable, intent(out) :: x_flux(:, :, :), y_flux(:, :, :)
    integer, intent(inout) :: level_advances(:)
    character(len=*), intent(inout) :: failure_context
    logical, intent(out) :: ok
    type(reactive_eb_exterior_state_2d), intent(in), optional :: exterior

    type(eb_geometry_2d) :: geometry, child_geometry
    type(amr_eb_flux_register_2d), allocatable :: registers(:)
    type(reactive_eb_exterior_state_2d) :: child_exterior
    real(dp), allocatable :: state_start(:, :, :), temperature_start(:, :)
    real(dp), allocatable :: state_end(:, :, :), temperature_end(:, :)
    real(dp), allocatable :: parent_work(:, :, :), parent_work_temperature(:, :)
    real(dp), allocatable :: child_work(:, :, :), child_work_temperature(:, :)
    real(dp), allocatable :: child_x_flux(:, :, :), child_y_flux(:, :, :)
    real(dp), allocatable :: integral_before(:)
    real(dp) :: alpha, child_interval
    integer :: child, child_count, first_child, global_child, ratio, substep
    logical :: local_ok, requires_closure

    ok = .false.
    if (.not. ieee_is_finite(interval) .or. interval <= 0.0_dp .or. &
        level < 1 .or. level > solution%level_count() .or. &
        patch_index < 1 .or. &
        patch_index > solution%levels(level)%patch_count() .or. &
        size(level_advances) /= solution%level_count()) return
    call patch_geometry_at( &
      solution%topology, level, patch_index, geometry, local_ok)
    if (.not. local_ok) return

    child_count = 0
    first_child = 0
    if (level < solution%level_count()) then
      first_child = solution%topology%relations(level)% &
        child_offsets(patch_index) + 1
      child_count = solution%topology%relations(level)% &
          child_offsets(patch_index + 1) - &
        solution%topology%relations(level)%child_offsets(patch_index)
    end if
    requires_closure = child_count > 0
    if (requires_closure) then
      allocate(integral_before(solution%nvar))
      call composite_reactive_amr_eb_patch_subtree_integral_2d( &
        solution, level, patch_index, integral_before, local_ok)
      if (.not. local_ok) return
    end if

    allocate(state_start, source= &
      solution%levels(level)%patches(patch_index)%state)
    allocate(temperature_start, source= &
      solution%levels(level)%patches(patch_index)%temperature)
    allocate(state_end, mold=state_start)
    allocate(temperature_end, mold=temperature_start)
    allocate(x_flux(solution%nvar, 0:geometry%nx, geometry%ny))
    allocate(y_flux(solution%nvar, geometry%nx, 0:geometry%ny))
    write(failure_context, '(a,i0,a,i0)') &
      "level advance level ", level - 1, " patch ", patch_index
    if (present(exterior)) then
      call advance_reactive_eb_level_2d( &
        species, state_start, temperature_start, geometry, solver, &
        reconstruction, limiter, selected_target, state_redist_max_order, &
        interval, state_end, temperature_end, x_flux, y_flux, local_ok, &
        exterior)
    else
      call advance_reactive_eb_level_2d( &
        species, state_start, temperature_start, geometry, solver, &
        reconstruction, limiter, selected_target, state_redist_max_order, &
        interval, state_end, temperature_end, x_flux, y_flux, local_ok)
    end if
    if (.not. local_ok) return
    solution%levels(level)%patches(patch_index)%state = state_end
    solution%levels(level)%patches(patch_index)%temperature = temperature_end
    level_advances(level) = level_advances(level) + 1
    if (child_count == 0) then
      ok = .true.
      return
    end if

    ratio = solution%topology%relations(level)%refinement_ratio
    child_interval = interval / real(ratio, dp)
    allocate(registers(child_count))
    do child = 1, child_count
      global_child = first_child + child - 1
      child_geometry = solution%topology%relations(level)% &
        children(global_child)%geometry
      write(failure_context, '(a,i0,a,i0)') &
        "flux register level ", level, " child ", global_child
      call initialize_amr_eb_flux_register_2d( &
        geometry, child_geometry, &
        solution%topology%relations(level)%children(global_child)%patch, &
        solution%nvar, registers(child), local_ok)
      if (.not. local_ok) return
      call accumulate_coarse_eb_fluxes_2d( &
        registers(child), geometry, child_geometry, &
        solution%topology%relations(level)%children(global_child)%patch, &
        x_flux, y_flux, interval, local_ok)
      if (.not. local_ok) return
    end do

    do substep = 1, ratio
      alpha = patch_tree_substep_time_alpha(reconstruction, substep, ratio)
      do child = 1, child_count
        global_child = first_child + child - 1
        child_geometry = solution%topology%relations(level)% &
          children(global_child)%geometry
        write(failure_context, '(a,i0,a,i0,a,i0)') &
          "exterior level ", level, " child ", global_child, &
          " substep ", substep
        call build_reactive_eb_patch_exterior_2d( &
          species, state_start, temperature_start, state_end, &
          temperature_end, geometry, child_geometry, &
          solution%topology%relations(level)%children(global_child)%patch, &
          alpha, child_exterior, local_ok, &
          solution%levels(level + 1)%patches(global_child)%state, &
          solution%levels(level + 1)%patches(global_child)%temperature)
        if (.not. local_ok) return
        call advance_reactive_amr_eb_patch_tree_hydro_node_2d( &
          species, solution, level + 1, global_child, solver, &
          reconstruction, limiter, state_redist_max_order, selected_target, &
          child_interval, child_x_flux, child_y_flux, level_advances, &
          failure_context, local_ok, child_exterior)
        if (.not. local_ok) return
        write(failure_context, '(a,i0,a,i0,a,i0)') &
          "fine flux level ", level, " child ", global_child, &
          " substep ", substep
        call accumulate_fine_eb_fluxes_2d( &
          registers(child), geometry, child_geometry, &
          solution%topology%relations(level)%children(global_child)%patch, &
          child_x_flux, child_y_flux, child_interval, local_ok)
        if (.not. local_ok) return
      end do
    end do

    allocate(parent_work, mold=state_end)
    allocate(parent_work_temperature, mold=temperature_end)
    do child = 1, child_count
      global_child = first_child + child - 1
      child_geometry = solution%topology%relations(level)% &
        children(global_child)%geometry
      if (allocated(child_work)) deallocate(child_work)
      if (allocated(child_work_temperature)) deallocate(child_work_temperature)
      allocate(child_work, mold= &
        solution%levels(level + 1)%patches(global_child)%state)
      allocate(child_work_temperature, mold= &
        solution%levels(level + 1)%patches(global_child)%temperature)
      write(failure_context, '(a,i0,a,i0)') &
        "reflux level ", level, " child ", global_child
      call reflux_reactive_eb_state_patch_2d( &
        species, solution%levels(level)%patches(patch_index)%state, &
        solution%levels(level)%patches(patch_index)%temperature, geometry, &
        solution%levels(level + 1)%patches(global_child)%state, &
        solution%levels(level + 1)%patches(global_child)%temperature, &
        child_geometry, &
        solution%topology%relations(level)%children(global_child)%patch, &
        registers(child), parent_work, parent_work_temperature, child_work, &
        child_work_temperature, local_ok)
      if (.not. local_ok) return
      solution%levels(level)%patches(patch_index)%state = parent_work
      solution%levels(level)%patches(patch_index)%temperature = &
        parent_work_temperature
      solution%levels(level + 1)%patches(global_child)%state = child_work
      solution%levels(level + 1)%patches(global_child)%temperature = &
        child_work_temperature
    end do

    do child = 1, child_count
      global_child = first_child + child - 1
      child_geometry = solution%topology%relations(level)% &
        children(global_child)%geometry
      write(failure_context, '(a,i0,a,i0)') &
        "average down level ", level, " child ", global_child
      call average_down_reactive_eb_state_patch_2d( &
        species, solution%levels(level)%patches(patch_index)%state, &
        solution%levels(level)%patches(patch_index)%temperature, geometry, &
        solution%levels(level + 1)%patches(global_child)%state, &
        child_geometry, &
        solution%topology%relations(level)%children(global_child)%patch, &
        parent_work, parent_work_temperature, local_ok)
      if (.not. local_ok) return
      solution%levels(level)%patches(patch_index)%state = parent_work
      solution%levels(level)%patches(patch_index)%temperature = &
        parent_work_temperature
    end do

    if (requires_closure) then
      write(failure_context, '(a,i0,a,i0)') &
        "cut-interface closure level ", level - 1, " patch ", patch_index
      call close_reactive_amr_eb_patch_subtree_conservation_2d( &
        species, solution, level, patch_index, integral_before, x_flux, &
        y_flux, interval, local_ok)
      if (.not. local_ok) return
    end if
    ok = .true.
  end subroutine advance_reactive_amr_eb_patch_tree_hydro_node_2d

  subroutine composite_integral_reactive_amr_eb_patch_tree_2d( &
      solution, integral, ok)
    type(reactive_amr_eb_patch_tree_2d), intent(in) :: solution
    real(dp), intent(out) :: integral(:)
    logical, intent(out) :: ok

    type(eb_geometry_2d) :: parent_geometry
    type(amr_eb_patch_2d) :: patch
    integer :: child, component, i, j, parent, relation

    integral = 0.0_dp
    ok = .false.
    if (.not. solution%is_valid()) return
    if (size(integral) /= solution%nvar) return

    do component = 1, solution%nvar
      integral(component) = sum( &
        solution%topology%root_geometry%volume_fraction * &
        solution%levels(1)%patches(1)%state(component, :, :)) * &
        solution%topology%root_geometry%dx * &
        solution%topology%root_geometry%dy
    end do

    do relation = 1, size(solution%topology%relations)
      do child = 1, &
          solution%topology%relations(relation)%child_patch_count()
        parent = solution%topology%relations(relation)% &
          children(child)%parent_patch
        call patch_geometry_at( &
          solution%topology, relation, parent, parent_geometry, ok)
        if (.not. ok) return
        patch = solution%topology%relations(relation)%children(child)%patch
        do component = 1, solution%nvar
          integral(component) = integral(component) + sum( &
            solution%topology%relations(relation)%children(child)% &
              geometry%volume_fraction * &
            solution%levels(relation + 1)%patches(child)% &
              state(component, :, :)) * &
            solution%topology%relations(relation)%children(child)% &
              geometry%dx * &
            solution%topology%relations(relation)%children(child)% &
              geometry%dy
          do j = patch%coarse_j_lower, patch%coarse_j_upper
            do i = patch%coarse_i_lower, patch%coarse_i_upper
              integral(component) = integral(component) - &
                parent_geometry%volume_fraction(i, j) * &
                solution%levels(relation)%patches(parent)% &
                  state(component, i, j) * &
                parent_geometry%dx * parent_geometry%dy
            end do
          end do
        end do
      end do
    end do
    ok = all(ieee_is_finite(integral))
  end subroutine composite_integral_reactive_amr_eb_patch_tree_2d

  subroutine synchronize_candidate(species, candidate, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_amr_eb_patch_tree_2d), intent(inout) :: candidate
    logical, intent(out) :: ok

    type(eb_geometry_2d) :: geometry, parent_geometry
    real(dp), allocatable :: state_work(:, :, :), temperature_work(:, :)
    logical :: local_ok
    integer :: child, level, parent, patch, relation

    ok = .false.
    if (.not. candidate%is_valid()) return
    do level = 1, candidate%level_count()
      do patch = 1, candidate%levels(level)%patch_count()
        call patch_geometry_at( &
          candidate%topology, level, patch, geometry, local_ok)
        if (.not. local_ok) return
        call recover_patch_temperature( &
          species, candidate%levels(level)%patches(patch), geometry, local_ok)
        if (.not. local_ok) return
      end do
    end do
    do relation = size(candidate%topology%relations), 1, -1
      do child = 1, &
          candidate%topology%relations(relation)%child_patch_count()
        parent = candidate%topology%relations(relation)% &
          children(child)%parent_patch
        call patch_geometry_at( &
          candidate%topology, relation, parent, parent_geometry, local_ok)
        if (.not. local_ok) return
        allocate(state_work, mold= &
          candidate%levels(relation)%patches(parent)%state)
        allocate(temperature_work, mold= &
          candidate%levels(relation)%patches(parent)%temperature)
        call average_down_reactive_eb_state_patch_2d( &
          species, candidate%levels(relation)%patches(parent)%state, &
          candidate%levels(relation)%patches(parent)%temperature, &
          parent_geometry, &
          candidate%levels(relation + 1)%patches(child)%state, &
          candidate%topology%relations(relation)%children(child)%geometry, &
          candidate%topology%relations(relation)%children(child)%patch, &
          state_work, temperature_work, local_ok)
        if (.not. local_ok) return
        candidate%levels(relation)%patches(parent)%state = state_work
        candidate%levels(relation)%patches(parent)%temperature = &
          temperature_work
        deallocate(state_work, temperature_work)
      end do
    end do
    ok = candidate%is_valid()
  end subroutine synchronize_candidate

  subroutine retain_same_resolution_overlap( &
      new_node, new_geometry, old_node, old_geometry, copied, ok)
    type(reactive_amr_eb_patch_tree_node_2d), intent(inout) :: new_node
    type(eb_geometry_2d), intent(in) :: new_geometry
    type(reactive_amr_eb_patch_tree_node_2d), intent(in) :: old_node
    type(eb_geometry_2d), intent(in) :: old_geometry
    logical, intent(inout) :: copied(:, :)
    logical, intent(out) :: ok

    real(dp) :: offset_i_real, offset_j_real, spacing_scale
    integer :: i, j, offset_i, offset_j, old_i, old_j

    ok = .true.
    spacing_scale = max(1.0_dp, abs(new_geometry%dx), &
      abs(new_geometry%dy), abs(old_geometry%dx), abs(old_geometry%dy))
    if (abs(new_geometry%dx - old_geometry%dx) > &
          geometry_tolerance * spacing_scale .or. &
        abs(new_geometry%dy - old_geometry%dy) > &
          geometry_tolerance * spacing_scale) return
    offset_i_real = (new_geometry%x_lower - old_geometry%x_lower) / &
      old_geometry%dx
    offset_j_real = (new_geometry%y_lower - old_geometry%y_lower) / &
      old_geometry%dy
    offset_i = nint(offset_i_real)
    offset_j = nint(offset_j_real)
    if (abs(offset_i_real - real(offset_i, dp)) > geometry_tolerance .or. &
        abs(offset_j_real - real(offset_j, dp)) > geometry_tolerance) return

    do j = 1, new_geometry%ny
      old_j = j + offset_j
      if (old_j < 1 .or. old_j > old_geometry%ny) cycle
      do i = 1, new_geometry%nx
        if (copied(i, j)) cycle
        old_i = i + offset_i
        if (old_i < 1 .or. old_i > old_geometry%nx) cycle
        if (.not. overlap_cell_geometry_matches( &
            new_geometry, i, j, old_geometry, old_i, old_j)) then
          ok = .false.
          return
        end if
        new_node%state(:, i, j) = old_node%state(:, old_i, old_j)
        new_node%temperature(i, j) = old_node%temperature(old_i, old_j)
        copied(i, j) = .true.
      end do
    end do
  end subroutine retain_same_resolution_overlap

  pure logical function overlap_cell_geometry_matches( &
      first, first_i, first_j, second, second_i, second_j) result(matches)
    type(eb_geometry_2d), intent(in) :: first, second
    integer, intent(in) :: first_i, first_j, second_i, second_j

    matches = first%cell_type(first_i, first_j) == &
        second%cell_type(second_i, second_j) .and. &
      abs(first%volume_fraction(first_i, first_j) - &
        second%volume_fraction(second_i, second_j)) <= &
        geometry_tolerance .and. &
      abs(first%cell_centroid_x(first_i, first_j) - &
        second%cell_centroid_x(second_i, second_j)) <= &
        geometry_tolerance .and. &
      abs(first%cell_centroid_y(first_i, first_j) - &
        second%cell_centroid_y(second_i, second_j)) <= &
        geometry_tolerance .and. &
      abs(first%boundary_length(first_i, first_j) - &
        second%boundary_length(second_i, second_j)) <= &
        geometry_tolerance .and. &
      abs(first%boundary_centroid_x(first_i, first_j) - &
        second%boundary_centroid_x(second_i, second_j)) <= &
        geometry_tolerance .and. &
      abs(first%boundary_centroid_y(first_i, first_j) - &
        second%boundary_centroid_y(second_i, second_j)) <= &
        geometry_tolerance .and. &
      abs(first%boundary_normal_x(first_i, first_j) - &
        second%boundary_normal_x(second_i, second_j)) <= &
        geometry_tolerance .and. &
      abs(first%boundary_normal_y(first_i, first_j) - &
        second%boundary_normal_y(second_i, second_j)) <= &
        geometry_tolerance .and. &
      abs(first%boundary_normal_integral_x(first_i, first_j) - &
        second%boundary_normal_integral_x(second_i, second_j)) <= &
        geometry_tolerance .and. &
      abs(first%boundary_normal_integral_y(first_i, first_j) - &
        second%boundary_normal_integral_y(second_i, second_j)) <= &
        geometry_tolerance .and. &
      all(abs(first%x_face_fraction(first_i - 1:first_i, first_j) - &
        second%x_face_fraction(second_i - 1:second_i, second_j)) <= &
        geometry_tolerance) .and. &
      all(abs(first%x_face_centroid_y(first_i - 1:first_i, first_j) - &
        second%x_face_centroid_y(second_i - 1:second_i, second_j)) <= &
        geometry_tolerance) .and. &
      all(abs(first%y_face_fraction(first_i, first_j - 1:first_j) - &
        second%y_face_fraction(second_i, second_j - 1:second_j)) <= &
        geometry_tolerance) .and. &
      all(abs(first%y_face_centroid_x(first_i, first_j - 1:first_j) - &
        second%y_face_centroid_x(second_i, second_j - 1:second_j)) <= &
        geometry_tolerance)
  end function overlap_cell_geometry_matches

  subroutine recover_patch_temperature(species, node, geometry, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_amr_eb_patch_tree_node_2d), intent(inout) :: node
    type(eb_geometry_2d), intent(in) :: geometry
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:)
    real(dp) :: recovered_temperature, sound_speed
    logical :: local_ok
    integer :: i, j

    ok = .false.
    if (any(.not. ieee_is_finite(node%state)) .or. &
        any(.not. ieee_is_finite(node%temperature)) .or. &
        any(node%temperature <= 0.0_dp)) return
    allocate(primitive(reactive_nprim(size(species))))
    do j = 1, geometry%ny
      do i = 1, geometry%nx
        if (geometry%cell_type(i, j) == eb_covered_cell) cycle
        call reactive_conserved_to_primitive( &
          species, node%state(:, i, j), node%temperature(i, j), &
          primitive, recovered_temperature, sound_speed, local_ok)
        if (.not. local_ok) return
        node%temperature(i, j) = recovered_temperature
      end do
    end do
    ok = .true.
  end subroutine recover_patch_temperature

  recursive subroutine composite_reactive_amr_eb_patch_subtree_integral_2d( &
      solution, level, patch_index, integral, ok)
    type(reactive_amr_eb_patch_tree_2d), intent(in) :: solution
    integer, intent(in) :: level, patch_index
    real(dp), intent(out) :: integral(:)
    logical, intent(out) :: ok

    type(eb_geometry_2d) :: geometry
    type(amr_eb_patch_2d) :: child_patch
    real(dp), allocatable :: child_integral(:)
    logical, allocatable :: refined(:, :)
    integer :: child, first_child, global_child, i, j, last_child
    logical :: local_ok

    integral = 0.0_dp
    ok = .false.
    if (size(integral) /= solution%nvar .or. &
        level < 1 .or. level > solution%level_count() .or. &
        patch_index < 1 .or. &
        patch_index > solution%levels(level)%patch_count()) return
    call patch_geometry_at( &
      solution%topology, level, patch_index, geometry, local_ok)
    if (.not. local_ok) return
    allocate(refined(geometry%nx, geometry%ny), source=.false.)

    first_child = 1
    last_child = 0
    if (level < solution%level_count()) then
      first_child = solution%topology%relations(level)% &
        child_offsets(patch_index) + 1
      last_child = solution%topology%relations(level)% &
        child_offsets(patch_index + 1)
      do global_child = first_child, last_child
        child_patch = solution%topology%relations(level)% &
          children(global_child)%patch
        refined(child_patch%coarse_i_lower:child_patch%coarse_i_upper, &
          child_patch%coarse_j_lower:child_patch%coarse_j_upper) = .true.
      end do
    end if

    do j = 1, geometry%ny
      do i = 1, geometry%nx
        if (refined(i, j)) cycle
        integral = integral + geometry%volume_fraction(i, j) * &
          solution%levels(level)%patches(patch_index)%state(:, i, j) * &
          geometry%dx * geometry%dy
      end do
    end do
    allocate(child_integral(solution%nvar))
    do global_child = first_child, last_child
      call composite_reactive_amr_eb_patch_subtree_integral_2d( &
        solution, level + 1, global_child, child_integral, local_ok)
      if (.not. local_ok) return
      integral = integral + child_integral
    end do
    if (any(.not. ieee_is_finite(integral))) then
      integral = 0.0_dp
      return
    end if
    ok = .true.
  end subroutine composite_reactive_amr_eb_patch_subtree_integral_2d

  subroutine close_reactive_amr_eb_patch_subtree_conservation_2d( &
      species, solution, level, patch_index, integral_before, x_flux, &
      y_flux, interval, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_amr_eb_patch_tree_2d), intent(inout) :: solution
    integer, intent(in) :: level, patch_index
    real(dp), intent(in) :: integral_before(:)
    real(dp), intent(in) :: x_flux(:, 0:, :), y_flux(:, :, 0:)
    real(dp), intent(in) :: interval
    logical, intent(out) :: ok

    type(eb_geometry_2d) :: geometry
    real(dp), allocatable :: current_integral(:), boundary_change(:)
    real(dp), allocatable :: residual(:), correction(:)
    real(dp) :: recipient_volume, scale, closure_tolerance, species_residual
    logical :: local_ok
    integer :: component, i, j, k

    ok = .false.
    if (size(integral_before) /= solution%nvar .or. &
        size(x_flux, 1) /= solution%nvar .or. &
        size(y_flux, 1) /= solution%nvar .or. &
        .not. ieee_is_finite(interval) .or. interval <= 0.0_dp .or. &
        any(.not. ieee_is_finite(integral_before)) .or. &
        any(.not. ieee_is_finite(x_flux)) .or. &
        any(.not. ieee_is_finite(y_flux))) return
    call patch_geometry_at( &
      solution%topology, level, patch_index, geometry, local_ok)
    if (.not. local_ok) return
    if (any(shape(x_flux) /= &
          [solution%nvar, geometry%nx + 1, geometry%ny]) .or. &
        any(shape(y_flux) /= &
          [solution%nvar, geometry%nx, geometry%ny + 1])) return

    allocate(current_integral(solution%nvar), boundary_change(solution%nvar))
    allocate(residual(solution%nvar), correction(solution%nvar))
    call composite_reactive_amr_eb_patch_subtree_integral_2d( &
      solution, level, patch_index, current_integral, local_ok)
    if (.not. local_ok) return
    boundary_change = 0.0_dp
    do j = 1, geometry%ny
      boundary_change = boundary_change + interval * geometry%dy * &
        (geometry%x_face_fraction(0, j) * x_flux(:, 0, j) - &
         geometry%x_face_fraction(geometry%nx, j) * &
           x_flux(:, geometry%nx, j))
    end do
    do i = 1, geometry%nx
      boundary_change = boundary_change + interval * geometry%dx * &
        (geometry%y_face_fraction(i, 0) * y_flux(:, i, 0) - &
         geometry%y_face_fraction(i, geometry%ny) * &
           y_flux(:, i, geometry%ny))
    end do
    residual = integral_before + boundary_change - current_integral
    correction = 0.0_dp
    correction(irho) = residual(irho)
    correction(iet) = residual(iet)
    species_residual = 0.0_dp
    do k = 1, size(species)
      component = reactive_species_component(k)
      correction(component) = residual(component)
      species_residual = species_residual + residual(component)
    end do
    scale = max(1.0_dp, maxval(abs(integral_before)), &
      maxval(abs(current_integral)), maxval(abs(boundary_change)))
    closure_tolerance = 4096.0_dp * epsilon(1.0_dp) * scale
    if (abs(residual(irho) - species_residual) > closure_tolerance) return
    component = reactive_species_component(size(species))
    correction(component) = correction(component) + &
      residual(irho) - species_residual
    if (maxval(abs(correction)) <= closure_tolerance) then
      ok = .true.
      return
    end if

    recipient_volume = 0.0_dp
    do j = 1, geometry%ny
      do i = 1, geometry%nx
        if (patch_tree_parent_cell_is_refined( &
              solution, level, patch_index, i, j) .or. &
            geometry%cell_type(i, j) == eb_covered_cell) cycle
        recipient_volume = recipient_volume + &
          geometry%volume_fraction(i, j) * geometry%dx * geometry%dy
      end do
    end do
    if (.not. ieee_is_finite(recipient_volume) .or. &
        recipient_volume <= tiny(1.0_dp)) return
    correction = correction / recipient_volume
    if (any(.not. ieee_is_finite(correction))) return
    do j = 1, geometry%ny
      do i = 1, geometry%nx
        if (patch_tree_parent_cell_is_refined( &
              solution, level, patch_index, i, j) .or. &
            geometry%cell_type(i, j) == eb_covered_cell) cycle
        solution%levels(level)%patches(patch_index)%state(:, i, j) = &
          solution%levels(level)%patches(patch_index)%state(:, i, j) + &
            correction
      end do
    end do
    call recover_patch_temperature( &
      species, solution%levels(level)%patches(patch_index), geometry, local_ok)
    if (.not. local_ok) return

    call composite_reactive_amr_eb_patch_subtree_integral_2d( &
      solution, level, patch_index, current_integral, local_ok)
    if (.not. local_ok) return
    residual = integral_before + boundary_change - current_integral
    if (abs(residual(irho)) > 8.0_dp * closure_tolerance .or. &
        abs(residual(iet)) > 8.0_dp * closure_tolerance) return
    do k = 1, size(species)
      if (abs(residual(reactive_species_component(k))) > &
          8.0_dp * closure_tolerance) return
    end do
    ok = .true.
  end subroutine close_reactive_amr_eb_patch_subtree_conservation_2d

  pure logical function patch_tree_parent_cell_is_refined( &
      solution, level, patch_index, i, j) result(refined)
    type(reactive_amr_eb_patch_tree_2d), intent(in) :: solution
    integer, intent(in) :: level, patch_index, i, j

    type(amr_eb_patch_2d) :: patch
    integer :: child, first_child, last_child

    refined = .false.
    if (level >= solution%level_count()) return
    first_child = solution%topology%relations(level)% &
      child_offsets(patch_index) + 1
    last_child = solution%topology%relations(level)% &
      child_offsets(patch_index + 1)
    do child = first_child, last_child
      patch = solution%topology%relations(level)%children(child)%patch
      refined = i >= patch%coarse_i_lower .and. &
        i <= patch%coarse_i_upper .and. j >= patch%coarse_j_lower .and. &
        j <= patch%coarse_j_upper
      if (refined) return
    end do
  end function patch_tree_parent_cell_is_refined

  pure real(dp) function patch_tree_substep_time_alpha( &
      reconstruction, substep, ratio) result(alpha)
    character(len=*), intent(in) :: reconstruction
    integer, intent(in) :: substep, ratio

    if (trim(reconstruction) == "characteristic_plm") then
      alpha = (real(substep, dp) - 0.5_dp) / real(ratio, dp)
    else
      alpha = real(substep - 1, dp) / real(ratio, dp)
    end if
  end function patch_tree_substep_time_alpha

  pure logical function patch_tree_checkpoint_metadata_is_valid( &
      time, steps, regrids, minimum_dt) result(valid)
    real(dp), intent(in) :: time, minimum_dt
    integer, intent(in) :: steps, regrids

    real(dp) :: tolerance

    valid = ieee_is_finite(time) .and. ieee_is_finite(minimum_dt) .and. &
      time >= 0.0_dp .and. minimum_dt >= 0.0_dp .and. &
      steps >= 0 .and. regrids >= 0
    if (.not. valid) return
    if (steps == 0) then
      valid = time == 0.0_dp .and. minimum_dt == 0.0_dp .and. &
        regrids <= 1
      return
    end if
    if (steps == huge(steps)) then
      valid = .false.
      return
    end if
    tolerance = 64.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(time))
    valid = time > 0.0_dp .and. minimum_dt > 0.0_dp .and. &
      minimum_dt <= time + tolerance .and. regrids <= steps + 1
  end function patch_tree_checkpoint_metadata_is_valid

  subroutine write_patch_tree_checkpoint_geometry_2d( &
      unit, geometry, status)
    integer, intent(in) :: unit
    type(eb_geometry_2d), intent(in) :: geometry
    integer, intent(out) :: status

    status = 1
    if (.not. geometry%is_valid()) return
    if (geometry%ny < 1) return
    if (geometry%nx > checkpoint_maximum_geometry_cells / geometry%ny) return
    write(unit, '(2(i0,1x),6(es27.18e3,1x))', iostat=status) &
      geometry%nx, geometry%ny, geometry%x_lower, geometry%x_upper, &
      geometry%y_lower, geometry%y_upper, geometry%dx, geometry%dy
    if (status /= 0) return
    call write_patch_tree_checkpoint_real_field_2d( &
      unit, geometry%volume_fraction, status)
    call write_patch_tree_checkpoint_real_field_2d( &
      unit, geometry%cell_centroid_x, status)
    call write_patch_tree_checkpoint_real_field_2d( &
      unit, geometry%cell_centroid_y, status)
    call write_patch_tree_checkpoint_real_field_2d( &
      unit, geometry%x_face_fraction, status)
    call write_patch_tree_checkpoint_real_field_2d( &
      unit, geometry%y_face_fraction, status)
    call write_patch_tree_checkpoint_real_field_2d( &
      unit, geometry%x_face_centroid_y, status)
    call write_patch_tree_checkpoint_real_field_2d( &
      unit, geometry%y_face_centroid_x, status)
    call write_patch_tree_checkpoint_real_field_2d( &
      unit, geometry%boundary_length, status)
    call write_patch_tree_checkpoint_real_field_2d( &
      unit, geometry%boundary_centroid_x, status)
    call write_patch_tree_checkpoint_real_field_2d( &
      unit, geometry%boundary_centroid_y, status)
    call write_patch_tree_checkpoint_real_field_2d( &
      unit, geometry%boundary_normal_x, status)
    call write_patch_tree_checkpoint_real_field_2d( &
      unit, geometry%boundary_normal_y, status)
    call write_patch_tree_checkpoint_real_field_2d( &
      unit, geometry%boundary_normal_integral_x, status)
    call write_patch_tree_checkpoint_real_field_2d( &
      unit, geometry%boundary_normal_integral_y, status)
    call write_patch_tree_checkpoint_integer_field_2d( &
      unit, geometry%cell_type, status)
  end subroutine write_patch_tree_checkpoint_geometry_2d

  subroutine read_patch_tree_checkpoint_geometry_2d( &
      unit, geometry, status)
    integer, intent(in) :: unit
    type(eb_geometry_2d), intent(out) :: geometry
    integer, intent(out) :: status

    real(dp) :: bounds_and_spacing(6)
    integer :: nx, ny

    geometry = eb_geometry_2d()
    read(unit, *, iostat=status) nx, ny, bounds_and_spacing
    if (status /= 0) return
    if (nx < 1 .or. ny < 1) then
      status = 1
      return
    end if
    if (nx > checkpoint_maximum_geometry_cells / ny .or. &
        any(.not. ieee_is_finite(bounds_and_spacing))) then
      status = 1
      return
    end if

    geometry%nx = nx
    geometry%ny = ny
    geometry%x_lower = bounds_and_spacing(1)
    geometry%x_upper = bounds_and_spacing(2)
    geometry%y_lower = bounds_and_spacing(3)
    geometry%y_upper = bounds_and_spacing(4)
    geometry%dx = bounds_and_spacing(5)
    geometry%dy = bounds_and_spacing(6)
    allocate(geometry%volume_fraction(1:nx, 1:ny))
    allocate(geometry%cell_centroid_x(1:nx, 1:ny))
    allocate(geometry%cell_centroid_y(1:nx, 1:ny))
    allocate(geometry%x_face_fraction(0:nx, 1:ny))
    allocate(geometry%y_face_fraction(1:nx, 0:ny))
    allocate(geometry%x_face_centroid_y(0:nx, 1:ny))
    allocate(geometry%y_face_centroid_x(1:nx, 0:ny))
    allocate(geometry%boundary_length(1:nx, 1:ny))
    allocate(geometry%boundary_centroid_x(1:nx, 1:ny))
    allocate(geometry%boundary_centroid_y(1:nx, 1:ny))
    allocate(geometry%boundary_normal_x(1:nx, 1:ny))
    allocate(geometry%boundary_normal_y(1:nx, 1:ny))
    allocate(geometry%boundary_normal_integral_x(1:nx, 1:ny))
    allocate(geometry%boundary_normal_integral_y(1:nx, 1:ny))
    allocate(geometry%cell_type(1:nx, 1:ny))

    call read_patch_tree_checkpoint_real_field_2d( &
      unit, geometry%volume_fraction, status)
    call read_patch_tree_checkpoint_real_field_2d( &
      unit, geometry%cell_centroid_x, status)
    call read_patch_tree_checkpoint_real_field_2d( &
      unit, geometry%cell_centroid_y, status)
    call read_patch_tree_checkpoint_real_field_2d( &
      unit, geometry%x_face_fraction, status)
    call read_patch_tree_checkpoint_real_field_2d( &
      unit, geometry%y_face_fraction, status)
    call read_patch_tree_checkpoint_real_field_2d( &
      unit, geometry%x_face_centroid_y, status)
    call read_patch_tree_checkpoint_real_field_2d( &
      unit, geometry%y_face_centroid_x, status)
    call read_patch_tree_checkpoint_real_field_2d( &
      unit, geometry%boundary_length, status)
    call read_patch_tree_checkpoint_real_field_2d( &
      unit, geometry%boundary_centroid_x, status)
    call read_patch_tree_checkpoint_real_field_2d( &
      unit, geometry%boundary_centroid_y, status)
    call read_patch_tree_checkpoint_real_field_2d( &
      unit, geometry%boundary_normal_x, status)
    call read_patch_tree_checkpoint_real_field_2d( &
      unit, geometry%boundary_normal_y, status)
    call read_patch_tree_checkpoint_real_field_2d( &
      unit, geometry%boundary_normal_integral_x, status)
    call read_patch_tree_checkpoint_real_field_2d( &
      unit, geometry%boundary_normal_integral_y, status)
    call read_patch_tree_checkpoint_integer_field_2d( &
      unit, geometry%cell_type, status)
    if (status /= 0 .or. .not. geometry%is_valid()) status = 1
  end subroutine read_patch_tree_checkpoint_geometry_2d

  subroutine write_patch_tree_checkpoint_real_field_2d( &
      unit, field, status)
    integer, intent(in) :: unit
    real(dp), intent(in) :: field(:, :)
    integer, intent(inout) :: status

    integer :: j

    if (status /= 0) return
    do j = 1, size(field, 2)
      write(unit, '(*(es27.18e3,1x))', iostat=status) field(:, j)
      if (status /= 0) return
    end do
  end subroutine write_patch_tree_checkpoint_real_field_2d

  subroutine read_patch_tree_checkpoint_real_field_2d( &
      unit, field, status)
    integer, intent(in) :: unit
    real(dp), intent(out) :: field(:, :)
    integer, intent(inout) :: status

    integer :: j

    if (status /= 0) return
    do j = 1, size(field, 2)
      read(unit, *, iostat=status) field(:, j)
      if (status /= 0) return
    end do
  end subroutine read_patch_tree_checkpoint_real_field_2d

  subroutine write_patch_tree_checkpoint_integer_field_2d( &
      unit, field, status)
    integer, intent(in) :: unit
    integer, intent(in) :: field(:, :)
    integer, intent(inout) :: status

    integer :: j

    if (status /= 0) return
    do j = 1, size(field, 2)
      write(unit, '(*(i0,1x))', iostat=status) field(:, j)
      if (status /= 0) return
    end do
  end subroutine write_patch_tree_checkpoint_integer_field_2d

  subroutine read_patch_tree_checkpoint_integer_field_2d( &
      unit, field, status)
    integer, intent(in) :: unit
    integer, intent(out) :: field(:, :)
    integer, intent(inout) :: status

    integer :: j

    if (status /= 0) return
    do j = 1, size(field, 2)
      read(unit, *, iostat=status) field(:, j)
      if (status /= 0) return
    end do
  end subroutine read_patch_tree_checkpoint_integer_field_2d

  subroutine patch_geometry_at( &
      topology, level_index, patch_index, geometry, ok)
    type(amr_eb_patch_tree_topology_2d), intent(in) :: topology
    integer, intent(in) :: level_index, patch_index
    type(eb_geometry_2d), intent(out) :: geometry
    logical, intent(out) :: ok

    ok = .false.
    if (level_index < 1 .or. &
        level_index > topology%level_count()) return
    if (patch_index < 1 .or. &
        patch_index > topology%level_patch_count(level_index - 1)) return
    if (level_index == 1) then
      geometry = topology%root_geometry
    else
      geometry = topology%relations(level_index - 1)% &
        children(patch_index)%geometry
    end if
    ok = geometry%is_valid()
  end subroutine patch_geometry_at

end module amr_eb_patch_tree_reactive_2d_mod
