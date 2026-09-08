program test_mpi_amr_sparse_integrator_policy
  use, intrinsic :: iso_fortran_env, only: error_unit
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: &
    load_h2o2_elementary_mechanism
  use simulation_config_reactive_1d_mod, only: reactive_1d_config
  use amr_patch_tree_1d_mod, only: amr_patch_level_plan_1d
  use amr_patch_tree_reactive_1d_mod, only: &
    amr_patch_tree_reactive_solution_1d, initialize_patch_tree_reactive_1d
  use mpi_amr_patch_1d_mod, only: &
    mpi_amr_patch_distribution_1d, &
    initialize_mpi_amr_patch_distribution_1d
  use mpi_amr_sparse_patch_1d_mod, only: &
    mpi_amr_sparse_reactive_solution_1d, &
    scatter_owned_patch_tree_reactive_1d, &
    advance_sparse_patch_tree_chemistry_1d
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(amr_patch_level_plan_1d), allocatable :: empty_plans(:)
  type(amr_patch_tree_reactive_solution_1d) :: root_solution
  type(mpi_amr_patch_distribution_1d) :: distribution
  type(mpi_amr_sparse_reactive_solution_1d) :: solution, backup
  type(reactive_1d_config) :: config
  logical :: ok
  integer :: ierr, rank

  call MPI_Init(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Init failed"
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  if (ierr /= MPI_SUCCESS) call fail("MPI_Comm_rank failed")

  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok) call fail("MPI AMR integrator thermodynamics load failed")
  call load_h2o2_elementary_mechanism(reactions, ok)
  if (.not. ok) call fail("MPI AMR integrator mechanism load failed")
  call configure_case(config)
  allocate(empty_plans(0))
  call initialize_patch_tree_reactive_1d( &
    species, config, empty_plans, root_solution, ok)
  if (.not. ok) call fail("MPI AMR integrator root initialization failed")
  call initialize_mpi_amr_patch_distribution_1d( &
    root_solution%hierarchy, MPI_COMM_WORLD, distribution, ok)
  if (.not. ok) call fail("MPI AMR integrator distribution failed")
  call scatter_owned_patch_tree_reactive_1d( &
    distribution, root_solution, solution, ok)
  if (.not. ok) call fail("MPI AMR integrator scatter failed")
  backup = solution

  call advance_sparse_patch_tree_chemistry_1d( &
    species, reactions, config, 1.0e-12_dp, distribution, solution, ok, &
    chemistry_integrator="unsupported")
  if (ok .or. .not. same_owned_state(solution, backup)) then
    call fail("Invalid MPI AMR integrator did not roll back collectively")
  end if
  call advance_sparse_patch_tree_chemistry_1d( &
    species, reactions, config, 1.0e-12_dp, distribution, solution, ok, &
    chemistry_integrator="explicit")
  if (.not. ok .or. .not. solution%is_valid(distribution)) then
    call fail("Explicit MPI AMR integrator was not forwarded")
  end if

  if (rank == 0) &
    write(*, '(a)') "test_mpi_amr_sparse_integrator_policy: PASS"
  call MPI_Finalize(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Finalize failed"

contains

  subroutine configure_case(local_config)
    type(reactive_1d_config), intent(out) :: local_config

    local_config = reactive_1d_config()
    local_config%nx = 8
    local_config%x_lower = 0.0_dp
    local_config%x_upper = 0.008_dp
    local_config%final_time = 1.0e-9_dp
    local_config%problem = "uniform_reactor"
    local_config%boundary_condition = "periodic"
    local_config%chemistry_model = "elementary"
    local_config%chemistry_enabled = .true.
    local_config%transport_enabled = .false.
    local_config%initial_temperature = 1200.0_dp
    local_config%initial_pressure = 101325.0_dp
    local_config%x_h2 = 0.29570_dp
    local_config%x_h = 1.0e-5_dp
    local_config%x_o = 1.0e-5_dp
    local_config%x_o2 = 0.14784_dp
    local_config%x_oh = 1.0e-5_dp
    local_config%x_h2o = 0.0_dp
    local_config%x_n2 = 0.55643_dp
    local_config%amr_enabled = .true.
    local_config%amr_max_levels = 2
    local_config%amr_minimum_patch_cells = 2
  end subroutine configure_case

  pure logical function same_owned_state(left, right) result(same)
    type(mpi_amr_sparse_reactive_solution_1d), intent(in) :: left, right

    integer :: level, patch

    same = left%time == right%time .and. left%steps == right%steps .and. &
      all(left%level_advances == right%level_advances) .and. &
      all(left%transport_level_advances == right%transport_level_advances)
    if (.not. same) return
    do level = 1, size(left%levels)
      do patch = 1, size(left%levels(level)%patches)
        if (.not. left%levels(level)%is_local(patch)) cycle
        same = all(abs( &
          left%levels(level)%patches(patch)%state - &
          right%levels(level)%patches(patch)%state) <= 0.0_dp)
        if (.not. same) return
      end do
    end do
  end function same_owned_state

  subroutine fail(message)
    character(len=*), intent(in) :: message

    integer :: abort_ierr

    if (rank == 0) write(error_unit, '(a)') trim(message)
    call MPI_Abort(MPI_COMM_WORLD, 1, abort_ierr)
    error stop 1
  end subroutine fail

end program test_mpi_amr_sparse_integrator_policy
