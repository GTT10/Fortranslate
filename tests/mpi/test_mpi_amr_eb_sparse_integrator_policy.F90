program test_mpi_amr_eb_sparse_integrator_policy
  use, intrinsic :: iso_fortran_env, only: error_unit
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: &
    load_h2o2_elementary_mechanism
  use simulation_config_reactive_2d_mod, only: reactive_2d_config
  use reactive_1d_mod, only: reactive_nvar
  use reactive_2d_mod, only: initialize_reactive_2d
  use eb_geometry_2d_mod, only: &
    eb_geometry_2d, eb_covered_cell, eb_cut_cell, build_eb_geometry_2d
  use amr_eb_patch_tree_2d_mod, only: &
    amr_eb_patch_tree_level_plan_2d, amr_eb_patch_tree_topology_2d, &
    initialize_amr_eb_patch_tree_topology_2d
  use mpi_amr_eb_patch_tree_2d_mod, only: &
    mpi_amr_eb_patch_tree_distribution_2d, &
    mpi_sparse_reactive_amr_eb_patch_tree_2d, &
    initialize_mpi_amr_eb_patch_tree_distribution_2d, &
    initialize_sparse_owned_reactive_amr_eb_patch_tree_root_2d, &
    advance_sparse_owned_reactive_amr_eb_patch_tree_chemistry_2d
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(amr_eb_patch_tree_level_plan_2d), allocatable :: empty_plans(:)
  type(amr_eb_patch_tree_topology_2d) :: topology
  type(mpi_amr_eb_patch_tree_distribution_2d) :: distribution
  type(mpi_sparse_reactive_amr_eb_patch_tree_2d) :: solution, backup
  type(eb_geometry_2d) :: geometry
  type(reactive_2d_config) :: config
  real(dp), allocatable :: level_set(:, :)
  real(dp), allocatable :: root_state(:, :, :), root_temperature(:, :)
  real(dp) :: base_density, dx, dy
  character(len=32) :: rejected_policy
  logical :: ok
  integer :: ierr, i, j, nranks, rank

  call MPI_Init(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Init failed"
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  if (ierr /= MPI_SUCCESS) call fail("MPI_Comm_rank failed")
  call MPI_Comm_size(MPI_COMM_WORLD, nranks, ierr)
  if (ierr /= MPI_SUCCESS) call fail("MPI_Comm_size failed")

  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok) call fail("MPI EB integrator thermodynamics load failed")
  call load_h2o2_elementary_mechanism(reactions, ok)
  if (.not. ok) call fail("MPI EB integrator mechanism load failed")
  call configure_case(config)

  allocate(level_set(0:config%nx, 0:config%ny))
  do j = 0, config%ny
    do i = 0, config%nx
      level_set(i, j) = real(i, dp) - 2.5_dp
    end do
  end do
  call build_eb_geometry_2d( &
    level_set, config%x_lower, config%x_upper, config%y_lower, &
    config%y_upper, geometry, ok)
  if (.not. ok .or. .not. any(geometry%cell_type == eb_cut_cell) .or. &
      .not. any(geometry%cell_type == eb_covered_cell)) then
    call fail("MPI EB integrator geometry initialization failed")
  end if
  allocate(empty_plans(0))
  call initialize_amr_eb_patch_tree_topology_2d( &
    geometry, empty_plans, topology, ok)
  if (.not. ok) call fail("MPI EB integrator topology failed")
  call initialize_mpi_amr_eb_patch_tree_distribution_2d( &
    topology, MPI_COMM_WORLD, distribution, ok, 2)
  if (.not. ok) call fail("MPI EB integrator distribution failed")
  if (distribution%is_local(0, 1)) then
    call initialize_reactive_2d( &
      species, config, root_state, root_temperature, dx, dy, base_density, ok)
    if (.not. ok) call fail("MPI EB integrator root state failed")
  end if
  call initialize_sparse_owned_reactive_amr_eb_patch_tree_root_2d( &
    distribution, topology, reactive_nvar(size(species)), root_state, &
    root_temperature, solution, ok)
  if (.not. ok) call fail("MPI EB integrator sparse initialization failed")
  backup = solution

  if (nranks == 1) then
    rejected_policy = "unsupported"
  else if (rank == 0) then
    rejected_policy = "explicit"
  else
    rejected_policy = "implicit"
  end if
  call advance_sparse_owned_reactive_amr_eb_patch_tree_chemistry_2d( &
    species, reactions, distribution, solution, 1.0e-12_dp, 1.0e-6_dp, &
    1.0e-12_dp, ok, chemistry_integrator=trim(rejected_policy))
  if (ok .or. .not. same_owned_state(solution, backup)) then
    call fail("Rejected MPI EB integrator mutated the sparse state")
  end if

  call advance_sparse_owned_reactive_amr_eb_patch_tree_chemistry_2d( &
    species, reactions, distribution, solution, 1.0e-12_dp, 1.0e-6_dp, &
    1.0e-12_dp, ok, chemistry_integrator="explicit")
  if (.not. ok .or. .not. solution%is_valid(distribution)) then
    call fail("Explicit MPI EB integrator was not forwarded")
  end if

  if (rank == 0) &
    write(*, '(a)') "test_mpi_amr_eb_sparse_integrator_policy: PASS"
  call MPI_Finalize(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Finalize failed"

contains

  subroutine configure_case(local_config)
    type(reactive_2d_config), intent(out) :: local_config

    local_config = reactive_2d_config()
    local_config%nx = 6
    local_config%ny = 6
    local_config%x_lower = 0.0_dp
    local_config%x_upper = 0.006_dp
    local_config%y_lower = 0.0_dp
    local_config%y_upper = 0.006_dp
    local_config%problem = "uniform_reactor"
    local_config%chemistry_model = "elementary"
    local_config%chemistry_enabled = .true.
    local_config%transport_enabled = .false.
    local_config%initial_temperature = 1200.0_dp
    local_config%initial_pressure = 101325.0_dp
    local_config%initial_velocity_x = 0.0_dp
    local_config%initial_velocity_y = 0.0_dp
  end subroutine configure_case

  pure logical function same_owned_state(left, right) result(same)
    type(mpi_sparse_reactive_amr_eb_patch_tree_2d), intent(in) :: left, right

    integer :: level, patch

    same = left%nvar == right%nvar .and. &
      (allocated(left%levels) .eqv. allocated(right%levels))
    if (.not. same .or. .not. allocated(left%levels)) return
    same = size(left%levels) == size(right%levels)
    if (.not. same) return
    do level = 1, size(left%levels)
      same = allocated(left%levels(level)%patches) .eqv. &
        allocated(right%levels(level)%patches)
      if (.not. same) return
      if (.not. allocated(left%levels(level)%patches)) cycle
      same = size(left%levels(level)%patches) == &
        size(right%levels(level)%patches)
      if (.not. same) return
      do patch = 1, size(left%levels(level)%patches)
        same = allocated(left%levels(level)%patches(patch)%state) .eqv. &
          allocated(right%levels(level)%patches(patch)%state)
        if (.not. same) return
        same = allocated(left%levels(level)%patches(patch)%temperature) .eqv. &
          allocated(right%levels(level)%patches(patch)%temperature)
        if (.not. same) return
        if (allocated(left%levels(level)%patches(patch)%state)) then
          same = all(shape(left%levels(level)%patches(patch)%state) == &
            shape(right%levels(level)%patches(patch)%state))
          if (.not. same) return
          same = all(left%levels(level)%patches(patch)%state == &
            right%levels(level)%patches(patch)%state)
          if (.not. same) return
        end if
        if (allocated(left%levels(level)%patches(patch)%temperature)) then
          same = all(shape(left%levels(level)%patches(patch)%temperature) == &
            shape(right%levels(level)%patches(patch)%temperature))
          if (.not. same) return
          same = all(left%levels(level)%patches(patch)%temperature == &
            right%levels(level)%patches(patch)%temperature)
          if (.not. same) return
        end if
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

end program test_mpi_amr_eb_sparse_integrator_policy
