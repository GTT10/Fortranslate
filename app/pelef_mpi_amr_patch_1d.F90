program pelef_mpi_amr_patch_1d
  use mpi_f08
  use precision_mod, only: dp
  use state_indices_mod, only: irho
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: &
    load_h2o2_elementary_mechanism
  use simulation_config_reactive_1d_mod, only: reactive_1d_config
  use amr_patch_tree_1d_mod, only: &
    amr_patch_level_plan_1d, amr_patch_tree_hierarchy_1d, &
    amr_patch_tree_level_fields_1d, initialize_patch_tree_1d, &
    prolong_patch_tree_1d
  use amr_patch_tree_reactive_1d_mod, only: &
    amr_patch_tree_reactive_solution_1d, &
    initialize_patch_tree_reactive_1d, advance_patch_tree_chemistry, &
    patch_tree_reactive_timestep_1d, advance_patch_tree_reactive_hydro_1d, &
    patch_tree_reactive_integrals_1d
  use mpi_amr_patch_1d_mod, only: &
    mpi_amr_patch_distribution_1d, mpi_amr_level_halos_1d, &
    initialize_mpi_amr_patch_distribution_1d, &
    synchronize_owned_patch_tree_fields_1d, &
    exchange_owned_adjacent_patch_halos_1d, &
    synchronize_owned_patch_tree_reactive_1d, &
    advance_owned_patch_tree_chemistry_1d, &
    advance_owned_patch_tree_hydro_1d
  implicit none

  integer, parameter :: variable_count = 3
  integer, parameter :: halo_width = 4
  real(dp), parameter :: stale_value = -huge(1.0_dp)

  type(amr_patch_tree_hierarchy_1d) :: hierarchy, comparison_hierarchy
  type(amr_patch_level_plan_1d), allocatable :: plans(:)
  type(amr_patch_tree_level_fields_1d), allocatable :: fields(:)
  type(mpi_amr_level_halos_1d), allocatable :: halos(:)
  type(mpi_amr_patch_distribution_1d) :: distribution
  type(mpi_amr_patch_distribution_1d) :: comparison_distribution
  type(mpi_amr_patch_distribution_1d) :: reactive_distribution
  type(mpi_amr_patch_distribution_1d) :: adjacent_distribution
  type(amr_patch_tree_reactive_solution_1d) :: initial_reactive
  type(amr_patch_tree_reactive_solution_1d) :: serial_reactive
  type(amr_patch_tree_reactive_solution_1d) :: distributed_reactive
  type(amr_patch_tree_reactive_solution_1d) :: rejected_reactive
  type(amr_patch_tree_reactive_solution_1d) :: rejected_backup
  type(amr_patch_tree_reactive_solution_1d) :: serial_hydro
  type(amr_patch_tree_reactive_solution_1d) :: distributed_hydro
  type(amr_patch_tree_reactive_solution_1d) :: adjacent_initial
  type(amr_patch_tree_reactive_solution_1d) :: adjacent_serial
  type(amr_patch_tree_reactive_solution_1d) :: adjacent_distributed
  type(amr_patch_level_plan_1d), allocatable :: reactive_plans(:)
  type(amr_patch_level_plan_1d), allocatable :: adjacent_reactive_plans(:)
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(reactive_1d_config) :: reactive_config
  type(reactive_1d_config) :: adjacent_reactive_config
  real(dp) :: root(variable_count, 64)
  real(dp), allocatable :: initial_integral(:), final_integral(:)
  real(dp) :: reactive_difference, conservation_error, hydro_dt, adjacent_dt
  logical :: ok
  integer :: ierr, rank, nranks, level, patch, variable, cell
  integer :: parent, child, left_patch, right_patch, layer, cross_rank_faces
  integer :: local_chemistry_advances, global_chemistry_advances
  integer :: expected_patch_advances, corrupt_owner
  integer :: local_hydro_advances, global_hydro_advances
  integer :: expected_hydro_advances, cross_owner_hydro_faces

  call MPI_Init(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Init failed"
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Comm_rank failed"
  call MPI_Comm_size(MPI_COMM_WORLD, nranks, ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Comm_size failed"

  call build_test_hierarchy(.false., plans, hierarchy, ok)
  call assert_all(ok, "valid AMR patch tree", rank)
  call initialize_mpi_amr_patch_distribution_1d( &
    hierarchy, MPI_COMM_WORLD, distribution, ok)
  call assert_all(ok, "deterministic MPI AMR patch distribution", rank)
  call assert_all(distribution%rank == rank, "local rank metadata", rank)
  call assert_all(distribution%nranks == nranks, "rank-count metadata", rank)
  call assert_all( &
    sum(distribution%rank_patch_counts) == 9, &
    "every root/fine patch has exactly one owner", rank)
  call assert_all( &
    sum(distribution%rank_cell_counts) == 152, &
    "owned cell work is conserved", rank)

  cross_rank_faces = 0
  do parent = 1, hierarchy%relations(1)%parent_patch_count()
    do child = 1, hierarchy%relations(1)%child_sets(parent)%patch_count() - 1
      if (hierarchy%relations(1)%child_sets(parent)%patches(child)%fine%upper + &
          1 /= hierarchy%relations(1)%child_sets(parent)% &
          patches(child + 1)%fine%lower) cycle
      left_patch = hierarchy%relations(1)%child_index(parent, child)
      right_patch = hierarchy%relations(1)%child_index(parent, child + 1)
      if (distribution%owner_of(1, left_patch) /= &
          distribution%owner_of(1, right_patch)) &
        cross_rank_faces = cross_rank_faces + 1
    end do
  end do
  if (nranks > 1) then
    call assert_all(cross_rank_faces >= 1, &
      "at least one adjacent face crosses ranks", rank)
  end if

  root = 0.0_dp
  call prolong_patch_tree_1d(root, hierarchy, fields, ok)
  call assert_all(ok, "patch-field allocation", rank)
  do level = 0, hierarchy%level_count() - 1
    do patch = 1, hierarchy%level_patch_count(level)
      fields(level + 1)%patches(patch)%values = stale_value
      if (.not. distribution%is_local(level, patch)) cycle
      do cell = 1, size(fields(level + 1)%patches(patch)%values, 2)
        do variable = 1, variable_count
          fields(level + 1)%patches(patch)%values(variable, cell) = &
            expected_value(level, patch, variable, cell)
        end do
      end do
    end do
  end do

  call exchange_owned_adjacent_patch_halos_1d( &
    distribution, hierarchy, fields, halo_width, halos, ok)
  call assert_all(ok, "owner-authoritative adjacent halo exchange", rank)
  do parent = 1, hierarchy%relations(1)%parent_patch_count()
    do child = 1, hierarchy%relations(1)%child_sets(parent)%patch_count() - 1
      if (hierarchy%relations(1)%child_sets(parent)%patches(child)%fine%upper + &
          1 /= hierarchy%relations(1)%child_sets(parent)% &
          patches(child + 1)%fine%lower) cycle
      left_patch = hierarchy%relations(1)%child_index(parent, child)
      right_patch = hierarchy%relations(1)%child_index(parent, child + 1)
      call assert_all(halos(2)%patches(left_patch)%has_right, &
        "left sibling receives a right halo", rank)
      call assert_all(halos(2)%patches(right_patch)%has_left, &
        "right sibling receives a left halo", rank)
      do layer = 1, halo_width
        do variable = 1, variable_count
          call assert_all( &
            halos(2)%patches(left_patch)%right(variable, layer) == &
              expected_value(1, right_patch, variable, layer), &
            "right-source halo value", rank)
          cell = size(fields(2)%patches(left_patch)%values, 2) - layer + 1
          call assert_all( &
            halos(2)%patches(right_patch)%left(variable, layer) == &
              expected_value(1, left_patch, variable, cell), &
            "left-source halo value", rank)
        end do
      end do
    end do
  end do

  call synchronize_owned_patch_tree_fields_1d( &
    distribution, hierarchy, fields, ok)
  call assert_all(ok, "owner-authoritative full-patch synchronization", rank)
  do level = 0, hierarchy%level_count() - 1
    do patch = 1, hierarchy%level_patch_count(level)
      do cell = 1, size(fields(level + 1)%patches(patch)%values, 2)
        do variable = 1, variable_count
          call assert_all( &
            fields(level + 1)%patches(patch)%values(variable, cell) == &
              expected_value(level, patch, variable, cell), &
            "replicated patch equals owner state", rank)
        end do
      end do
    end do
  end do

  if (nranks > 1) then
    call build_test_hierarchy( &
      rank == nranks - 1, plans, comparison_hierarchy, ok)
    call assert_all(ok, "comparison AMR patch tree", rank)
    call initialize_mpi_amr_patch_distribution_1d( &
      comparison_hierarchy, MPI_COMM_WORLD, comparison_distribution, ok)
    call assert_all(.not. ok, &
      "rank-inconsistent hierarchy is rejected collectively", rank)
  end if

  call load_h2o2_elementary_thermo(species, ok)
  call assert_all(ok, "reactive AMR thermodynamics", rank)
  call load_h2o2_elementary_mechanism(reactions, ok)
  call assert_all(ok, "reactive AMR chemistry mechanism", rank)
  call configure_reactive_case(reactive_config)
  call build_reactive_plans(reactive_plans)
  call initialize_patch_tree_reactive_1d( &
    species, reactive_config, reactive_plans, initial_reactive, ok)
  call assert_all(ok .and. initial_reactive%is_valid(), &
    "four-level reactive AMR initialization", rank)
  call initialize_mpi_amr_patch_distribution_1d( &
    initial_reactive%hierarchy, MPI_COMM_WORLD, reactive_distribution, ok)
  call assert_all(ok, "reactive AMR owner distribution", rank)

  serial_reactive = initial_reactive
  distributed_reactive = initial_reactive
  allocate(initial_integral( &
    size(initial_reactive%levels(1)%patches(1)%state, 1)))
  allocate(final_integral(size(initial_integral)))
  call patch_tree_reactive_integrals_1d( &
    initial_reactive, initial_integral, ok)
  call assert_all(ok, "initial reactive composite integral", rank)
  call advance_patch_tree_chemistry( &
    species, reactions, reactive_config, 1.0e-10_dp, serial_reactive, ok)
  call assert_all(ok .and. serial_reactive%is_valid(), &
    "serial patch-tree chemistry reference", rank)
  call advance_owned_patch_tree_chemistry_1d( &
    species, reactions, reactive_config, 1.0e-10_dp, &
    reactive_distribution, distributed_reactive, ok, &
    local_chemistry_advances)
  call assert_all(ok .and. distributed_reactive%is_valid(), &
    "owner-only distributed patch-tree chemistry", rank)
  call assert_all( &
    local_chemistry_advances == &
      reactive_distribution%rank_patch_counts(rank + 1), &
    "only locally owned patches execute chemistry", rank)
  call MPI_Allreduce( &
    local_chemistry_advances, global_chemistry_advances, 1, MPI_INTEGER, &
    MPI_SUM, MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS, &
    "owner-only chemistry execution reduction", rank)
  expected_patch_advances = sum(reactive_distribution%rank_patch_counts)
  call assert_all(global_chemistry_advances == expected_patch_advances, &
    "every reactive patch advances exactly once globally", rank)
  reactive_difference = reactive_solution_difference( &
    distributed_reactive, serial_reactive)
  call assert_all(reactive_difference <= 5.0e-13_dp, &
    "distributed chemistry matches serial patch tree", rank)
  call assert_all( &
    reactive_solution_difference(distributed_reactive, initial_reactive) > &
      100.0_dp * epsilon(1.0_dp), &
    "owner-only chemistry changes the reactive state", rank)
  call patch_tree_reactive_integrals_1d( &
    distributed_reactive, final_integral, ok)
  call assert_all(ok, "distributed reactive composite integral", rank)
  conservation_error = maxval(abs( &
    final_integral(1:5) - initial_integral(1:5)) / &
    max(1.0_dp, abs(initial_integral(1:5))))
  call assert_all(conservation_error <= 3.0e-10_dp, &
    "owner-only chemistry conserves mass momentum energy", rank)

  rejected_reactive = initial_reactive
  corrupt_owner = reactive_distribution%owner_of(3, 1)
  if (rank == corrupt_owner) &
    rejected_reactive%levels(4)%patches(1)%state(irho, 1) = -1.0_dp
  call synchronize_owned_patch_tree_reactive_1d( &
    reactive_distribution, rejected_reactive, ok)
  call assert_all(ok, "corrupt owner state synchronization", rank)
  rejected_backup = rejected_reactive
  call advance_owned_patch_tree_chemistry_1d( &
    species, reactions, reactive_config, 1.0e-10_dp, &
    reactive_distribution, rejected_reactive, ok, &
    local_chemistry_advances)
  call assert_all(.not. ok .and. local_chemistry_advances == 0, &
    "owner chemistry failure is rejected globally", rank)
  call assert_all( &
    reactive_solution_difference(rejected_reactive, rejected_backup) == &
      0.0_dp, "global chemistry rollback is exact", rank)

  serial_hydro = initial_reactive
  distributed_hydro = initial_reactive
  call patch_tree_reactive_timestep_1d( &
    species, reactive_config, initial_reactive, hydro_dt, ok)
  call assert_all(ok .and. hydro_dt > 0.0_dp, &
    "four-level reactive hydro timestep", rank)
  hydro_dt = min(0.10_dp * hydro_dt, 2.0e-8_dp)
  call advance_patch_tree_reactive_hydro_1d( &
    species, reactive_config, hydro_dt, serial_hydro, ok)
  call assert_all(ok .and. serial_hydro%is_valid(), &
    "serial four-level hydro reference", rank)
  call advance_owned_patch_tree_hydro_1d( &
    species, reactive_config, hydro_dt, reactive_distribution, &
    distributed_hydro, ok, local_hydro_advances)
  call assert_all(ok .and. distributed_hydro%is_valid(), &
    "owner-only four-level hydro", rank)
  expected_hydro_advances = expected_owned_hydro_advances( &
    reactive_distribution, initial_reactive%hierarchy, rank)
  call assert_all(local_hydro_advances == expected_hydro_advances, &
    "four-level hydro executes on owners only", rank)
  call MPI_Allreduce( &
    local_hydro_advances, global_hydro_advances, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    global_hydro_advances == 33, &
    "four-level hydro global subcycle count", rank)
  call assert_all(all(distributed_hydro%level_advances == [1, 4, 12, 16]), &
    "four-level distributed hydro level accounting", rank)
  reactive_difference = reactive_solution_difference( &
    distributed_hydro, serial_hydro)
  call assert_all(reactive_difference <= 5.0e-13_dp, &
    "distributed four-level hydro matches serial", rank)
  call patch_tree_reactive_integrals_1d( &
    distributed_hydro, final_integral, ok)
  call assert_all(ok, "distributed hydro composite integral", rank)
  conservation_error = maxval(abs( &
    final_integral - initial_integral) / &
    max(1.0_dp, abs(initial_integral)))
  call assert_all(conservation_error <= 3.0e-10_dp, &
    "owner-only four-level hydro conservation", rank)

  adjacent_reactive_config = reactive_config
  adjacent_reactive_config%problem = "entropy_wave"
  adjacent_reactive_config%amr_reconstruction = "ppm"
  adjacent_reactive_config%ppm_contact_steepening = .false.
  adjacent_reactive_config%ppm_shock_flattening = .false.
  adjacent_reactive_config%amr_hybrid_weno = .false.
  call build_adjacent_reactive_plans(adjacent_reactive_plans)
  call initialize_patch_tree_reactive_1d( &
    species, adjacent_reactive_config, adjacent_reactive_plans, &
    adjacent_initial, ok)
  call assert_all(ok .and. adjacent_initial%is_valid(), &
    "adjacent reactive PPM initialization", rank)
  call initialize_mpi_amr_patch_distribution_1d( &
    adjacent_initial%hierarchy, MPI_COMM_WORLD, adjacent_distribution, ok)
  call assert_all(ok, "adjacent reactive owner distribution", rank)
  cross_owner_hydro_faces = 0
  do child = 1, adjacent_initial%hierarchy%relations(1)% &
      child_sets(1)%patch_count() - 1
    left_patch = adjacent_initial%hierarchy%relations(1)%child_index(1, child)
    right_patch = adjacent_initial%hierarchy%relations(1)% &
      child_index(1, child + 1)
    if (adjacent_distribution%owner_of(1, left_patch) /= &
        adjacent_distribution%owner_of(1, right_patch)) &
      cross_owner_hydro_faces = cross_owner_hydro_faces + 1
  end do
  if (nranks > 1) call assert_all(cross_owner_hydro_faces >= 1, &
    "adjacent reactive face crosses MPI owners", rank)
  serial_hydro = adjacent_initial
  distributed_hydro = adjacent_initial
  call patch_tree_reactive_integrals_1d( &
    adjacent_initial, initial_integral, ok)
  call assert_all(ok, "adjacent initial composite integral", rank)
  call patch_tree_reactive_timestep_1d( &
    species, adjacent_reactive_config, adjacent_initial, adjacent_dt, ok)
  call assert_all(ok .and. adjacent_dt > 0.0_dp, &
    "adjacent reactive hydro timestep", rank)
  adjacent_dt = min(0.10_dp * adjacent_dt, 2.0e-8_dp)
  call advance_patch_tree_reactive_hydro_1d( &
    species, adjacent_reactive_config, adjacent_dt, serial_hydro, ok)
  call assert_all(ok, "serial adjacent PPM hydro reference", rank)
  call advance_owned_patch_tree_hydro_1d( &
    species, adjacent_reactive_config, adjacent_dt, adjacent_distribution, &
    distributed_hydro, ok, local_hydro_advances)
  call assert_all(ok .and. distributed_hydro%is_valid(), &
    "owner-only adjacent PPM hydro", rank)
  expected_hydro_advances = expected_owned_hydro_advances( &
    adjacent_distribution, adjacent_initial%hierarchy, rank)
  call assert_all(local_hydro_advances == expected_hydro_advances, &
    "adjacent PPM hydro executes on owners only", rank)
  call MPI_Allreduce( &
    local_hydro_advances, global_hydro_advances, 1, MPI_INTEGER, MPI_SUM, &
    MPI_COMM_WORLD, ierr)
  call assert_all(ierr == MPI_SUCCESS .and. &
    global_hydro_advances == 13, &
    "adjacent PPM global subcycle count", rank)
  call assert_all(all(distributed_hydro%level_advances == [1, 12]), &
    "adjacent PPM distributed level accounting", rank)
  reactive_difference = reactive_solution_difference( &
    distributed_hydro, serial_hydro)
  call assert_all(reactive_difference <= 5.0e-13_dp, &
    "cross-owner adjacent PPM matches serial", rank)
  call patch_tree_reactive_integrals_1d( &
    distributed_hydro, final_integral, ok)
  call assert_all(ok, "adjacent distributed composite integral", rank)
  conservation_error = maxval(abs( &
    final_integral - initial_integral) / &
    max(1.0_dp, abs(initial_integral)))
  call assert_all(conservation_error <= 3.0e-10_dp, &
    "cross-owner adjacent PPM conservation", rank)

  if (rank == 0) write(*, '(a,i0,a)') &
    "pelef_mpi_amr_patch_1d: PASS (", nranks, " ranks)"
  call MPI_Finalize(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Finalize failed"

contains

  subroutine build_test_hierarchy( &
      shift_last_upper, local_plans, local_hierarchy, local_ok)
    logical, intent(in) :: shift_last_upper
    type(amr_patch_level_plan_1d), allocatable, intent(out) :: local_plans(:)
    type(amr_patch_tree_hierarchy_1d), intent(out) :: local_hierarchy
    logical, intent(out) :: local_ok

    integer, parameter :: lowers(8) = [5, 9, 17, 25, 33, 37, 49, 53]
    integer :: uppers(8)
    integer :: entry

    uppers = [8, 12, 24, 28, 36, 44, 52, 60]
    if (shift_last_upper) uppers(8) = 59
    allocate(local_plans(1))
    local_plans(1)%refinement_ratio = 2
    allocate(local_plans(1)%patches(8))
    do entry = 1, 8
      local_plans(1)%patches(entry)%parent_patch = 1
      local_plans(1)%patches(entry)%lower = lowers(entry)
      local_plans(1)%patches(entry)%upper = uppers(entry)
    end do
    call initialize_patch_tree_1d( &
      64, 0.0_dp, 1.0_dp, local_plans, local_hierarchy, local_ok)
  end subroutine build_test_hierarchy

  subroutine configure_reactive_case(local_config)
    type(reactive_1d_config), intent(out) :: local_config

    local_config = reactive_1d_config()
    local_config%nx = 32
    local_config%x_lower = 0.0_dp
    local_config%x_upper = 0.012_dp
    local_config%cfl = 0.20_dp
    local_config%problem = "reactive_hotspot"
    local_config%riemann_solver = "rusanov"
    local_config%limiter = "mc"
    local_config%boundary_condition = "periodic"
    local_config%chemistry_enabled = .true.
    local_config%transport_enabled = .false.
    local_config%chemistry_relative_tolerance = 1.0e-8_dp
    local_config%chemistry_absolute_tolerance = 1.0e-14_dp
    local_config%initial_temperature = 1200.0_dp
    local_config%initial_pressure = 101325.0_dp
    local_config%initial_velocity = 0.0_dp
    local_config%hotspot_temperature_rise = 200.0_dp
    local_config%hotspot_center = 0.006_dp
    local_config%hotspot_width = 0.0012_dp
    local_config%amr_enabled = .true.
    local_config%amr_reconstruction = "pcm"
  end subroutine configure_reactive_case

  subroutine build_reactive_plans(local_plans)
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
  end subroutine build_reactive_plans

  subroutine build_adjacent_reactive_plans(local_plans)
    type(amr_patch_level_plan_1d), allocatable, intent(out) :: local_plans(:)

    integer :: entry

    allocate(local_plans(1))
    local_plans(1)%refinement_ratio = 2
    allocate(local_plans(1)%patches(6))
    do entry = 1, 6
      local_plans(1)%patches(entry)%parent_patch = 1
      local_plans(1)%patches(entry)%lower = 4 * entry
      local_plans(1)%patches(entry)%upper = 4 * entry + 3
    end do
  end subroutine build_adjacent_reactive_plans

  integer function expected_owned_hydro_advances( &
      local_distribution, local_hierarchy, local_rank) result(count)
    type(mpi_amr_patch_distribution_1d), intent(in) :: local_distribution
    type(amr_patch_tree_hierarchy_1d), intent(in) :: local_hierarchy
    integer, intent(in) :: local_rank

    integer :: local_level, local_patch, multiplier

    count = 0
    multiplier = 1
    do local_level = 0, local_hierarchy%level_count() - 1
      if (local_level > 0) multiplier = multiplier * &
        local_hierarchy%relations(local_level)%refinement_ratio
      do local_patch = 1, local_hierarchy%level_patch_count(local_level)
        if (local_distribution%owner_of(local_level, local_patch) == &
            local_rank) count = count + multiplier
      end do
    end do
  end function expected_owned_hydro_advances

  real(dp) function reactive_solution_difference(first, second) result(error)
    type(amr_patch_tree_reactive_solution_1d), intent(in) :: first, second

    integer :: local_level, local_patch

    error = huge(1.0_dp)
    if (first%level_count() /= second%level_count()) return
    if (size(first%level_advances) /= size(second%level_advances)) return
    if (size(first%transport_level_advances) /= &
        size(second%transport_level_advances)) return
    error = abs(first%time - second%time) / max(1.0_dp, abs(second%time))
    error = max(error, real(abs(first%steps - second%steps), dp))
    error = max(error, real(maxval(abs( &
      first%level_advances - second%level_advances)), dp))
    error = max(error, real(maxval(abs( &
      first%transport_level_advances - &
      second%transport_level_advances)), dp))
    do local_level = 1, first%level_count()
      if (size(first%levels(local_level)%patches) /= &
          size(second%levels(local_level)%patches)) then
        error = huge(1.0_dp)
        return
      end if
      do local_patch = 1, size(first%levels(local_level)%patches)
        error = max(error, maxval(abs( &
          first%levels(local_level)%patches(local_patch)%state - &
          second%levels(local_level)%patches(local_patch)%state) / &
          max(1.0_dp, abs(second%levels(local_level)% &
            patches(local_patch)%state))))
        error = max(error, maxval(abs( &
          first%levels(local_level)%patches(local_patch)%temperature - &
          second%levels(local_level)%patches(local_patch)%temperature) / &
          max(1.0_dp, abs(second%levels(local_level)% &
            patches(local_patch)%temperature))))
        error = max(error, maxval(abs( &
          first%levels(local_level)%patches(local_patch)%left_ghost_state - &
          second%levels(local_level)%patches(local_patch)%left_ghost_state) / &
          max(1.0_dp, abs(second%levels(local_level)% &
            patches(local_patch)%left_ghost_state))))
        error = max(error, maxval(abs( &
          first%levels(local_level)%patches(local_patch)%right_ghost_state - &
          second%levels(local_level)%patches(local_patch)% &
            right_ghost_state) / max(1.0_dp, abs(second%levels(local_level)% &
              patches(local_patch)%right_ghost_state))))
        error = max(error, maxval(abs( &
          first%levels(local_level)%patches(local_patch)% &
            left_ghost_temperature - second%levels(local_level)% &
            patches(local_patch)%left_ghost_temperature) / &
          max(1.0_dp, abs(second%levels(local_level)%patches(local_patch)% &
            left_ghost_temperature))))
        error = max(error, maxval(abs( &
          first%levels(local_level)%patches(local_patch)% &
            right_ghost_temperature - second%levels(local_level)% &
            patches(local_patch)%right_ghost_temperature) / &
          max(1.0_dp, abs(second%levels(local_level)%patches(local_patch)% &
            right_ghost_temperature))))
      end do
    end do
  end function reactive_solution_difference

  pure real(dp) function expected_value( &
      local_level, local_patch, local_variable, local_cell) result(value)
    integer, intent(in) :: local_level, local_patch, local_variable, local_cell

    value = real( &
      100000 * local_level + 10000 * local_patch + &
      100 * local_variable + local_cell, dp)
  end function expected_value

  subroutine assert_all(condition, message, local_rank)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    integer, intent(in) :: local_rank

    logical :: global_condition
    integer :: local_ierr

    call MPI_Allreduce(condition, global_condition, 1, MPI_LOGICAL, MPI_LAND, &
      MPI_COMM_WORLD, local_ierr)
    if (local_ierr == MPI_SUCCESS .and. global_condition) return
    if (local_rank == 0) write(*, '(a)') "FAIL: " // trim(message)
    call MPI_Abort(MPI_COMM_WORLD, 1, local_ierr)
  end subroutine assert_all

end program pelef_mpi_amr_patch_1d
