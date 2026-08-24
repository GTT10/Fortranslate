program pelef_mpi_amr_patch_1d
  use mpi_f08
  use precision_mod, only: dp
  use amr_patch_tree_1d_mod, only: &
    amr_patch_level_plan_1d, amr_patch_tree_hierarchy_1d, &
    amr_patch_tree_level_fields_1d, initialize_patch_tree_1d, &
    prolong_patch_tree_1d
  use mpi_amr_patch_1d_mod, only: &
    mpi_amr_patch_distribution_1d, mpi_amr_level_halos_1d, &
    initialize_mpi_amr_patch_distribution_1d, &
    synchronize_owned_patch_tree_fields_1d, &
    exchange_owned_adjacent_patch_halos_1d
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
  real(dp) :: root(variable_count, 64)
  logical :: ok
  integer :: ierr, rank, nranks, level, patch, variable, cell
  integer :: parent, child, left_patch, right_patch, layer, cross_rank_faces

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
