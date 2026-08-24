program pelef_mpi_amr_reactive_1d
  use, intrinsic :: iso_fortran_env, only: error_unit
  use mpi_f08
  use precision_mod, only: dp
  use constants_mod, only: pelef_version
  use state_indices_mod, only: irho, imx, imy, imz, iet
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use transport_database_mod, only: &
    gas_transport_species, load_h2o2_elementary_transport, &
    load_h2o2_full_transport
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use h2o2_elementary_mechanism_mod, only: &
    load_h2o2_elementary_mechanism
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use simulation_config_reactive_1d_mod, only: &
    reactive_1d_config, read_reactive_1d_configuration
  use reactive_1d_mod, only: reactive_nvar
  use amr_patch_tree_1d_mod, only: amr_patch_level_plan_1d
  use amr_patch_tree_reactive_1d_mod, only: &
    amr_patch_tree_reactive_solution_1d, &
    initialize_patch_tree_reactive_1d, patch_tree_reactive_integrals_1d, &
    write_patch_tree_reactive_1d_csv
  use mpi_amr_patch_1d_mod, only: &
    mpi_amr_patch_distribution_1d, &
    initialize_mpi_amr_patch_distribution_1d
  use mpi_amr_sparse_patch_1d_mod, only: &
    mpi_amr_sparse_reactive_solution_1d, &
    scatter_owned_patch_tree_reactive_1d, &
    materialize_owned_patch_tree_reactive_1d, &
    sparse_patch_tree_reactive_timestep_1d, &
    advance_sparse_patch_tree_reactive_1d, &
    regrid_tagged_sparse_patch_tree_reactive_1d
  implicit none

  type(reactive_1d_config) :: config
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(amr_patch_level_plan_1d), allocatable :: empty_plans(:)
  type(amr_patch_tree_reactive_solution_1d) :: root_solution
  type(amr_patch_tree_reactive_solution_1d) :: replicated_solution
  type(mpi_amr_sparse_reactive_solution_1d) :: sparse_solution
  type(mpi_amr_patch_distribution_1d) :: distribution, new_distribution
  real(dp), allocatable :: initial_all(:), final_all(:)
  real(dp) :: initial_integrals(5), final_integrals(5)
  real(dp) :: conservation_error(5), dt, tolerance
  character(len=1024) :: input_path, output_path, message
  logical :: ok, changed, output_ok
  integer :: ierr, rank, nranks, argument_count
  integer :: tagged_cells, transferred_cells

  call MPI_Init(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Init failed"
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Comm_rank failed"
  call MPI_Comm_size(MPI_COMM_WORLD, nranks, ierr)
  if (ierr /= MPI_SUCCESS) call abort_run("MPI_Comm_size failed", 2)

  argument_count = command_argument_count()
  if (argument_count < 1 .or. argument_count > 2) then
    call abort_run( &
      "Usage: pelef_mpi_amr_reactive_1d <input.nml> [output.csv]", 2)
  end if
  call get_command_argument(1, input_path)
  call read_reactive_1d_configuration(trim(input_path), config, ok, message)
  if (.not. ok) call abort_run(trim(message), 2)
  if (.not. config%amr_enabled) &
    call abort_run("Sparse MPI AMR requires amr_enabled", 2)
  if (config%amr_multipatch_enabled) &
    call abort_run( &
      "Sparse MPI AMR uses the patch-tree mode, not legacy multipatch mode", 2)
  output_path = config%output_file
  if (argument_count == 2) call get_command_argument(2, output_path)
  if (len_trim(output_path) == 0) &
    call abort_run("Sparse MPI AMR output path is empty", 2)

  select case (trim(config%chemistry_model))
  case ("elementary")
    call load_h2o2_elementary_thermo(species, ok)
    if (.not. ok) call abort_run("Failed to load elementary thermodynamics", 3)
    call load_h2o2_elementary_mechanism(reactions, ok)
    if (.not. ok) call abort_run("Failed to load elementary mechanism", 3)
    call load_h2o2_elementary_transport(transport, ok)
    if (.not. ok) call abort_run("Failed to load elementary transport", 3)
  case ("full_h2o2")
    call load_h2o2_full_thermo(species, ok)
    if (.not. ok) call abort_run("Failed to load full H2/O2 thermodynamics", 3)
    call load_h2o2_full_mechanism(reactions, ok)
    if (.not. ok) call abort_run("Failed to load full H2/O2 mechanism", 3)
    call load_h2o2_full_transport(transport, ok)
    if (.not. ok) call abort_run("Failed to load full H2/O2 transport", 3)
  case default
    call abort_run("Unknown chemistry model", 3)
  end select

  allocate(empty_plans(0))
  call initialize_patch_tree_reactive_1d( &
    species, config, empty_plans, root_solution, ok)
  if (.not. ok) call abort_run("Sparse MPI AMR root initialization failed", 4)
  allocate(initial_all(reactive_nvar(size(species))))
  allocate(final_all(reactive_nvar(size(species))))
  call patch_tree_reactive_integrals_1d(root_solution, initial_all, ok)
  if (.not. ok) call abort_run("Initial AMR integral failed", 4)
  initial_integrals = initial_all([irho, imx, imy, imz, iet])

  call initialize_mpi_amr_patch_distribution_1d( &
    root_solution%hierarchy, MPI_COMM_WORLD, distribution, ok, &
    config%amr_mpi_work_exponent)
  if (.not. ok) call abort_run("Initial sparse AMR ownership failed", 4)
  call scatter_owned_patch_tree_reactive_1d( &
    distribution, root_solution, sparse_solution, ok)
  if (.not. ok) call abort_run("Initial sparse AMR scatter failed", 4)
  call regrid_tagged_sparse_patch_tree_reactive_1d( &
    species, config, distribution, sparse_solution, new_distribution, &
    changed, tagged_cells, transferred_cells, ok)
  if (.not. ok) call abort_run("Initial sparse AMR tagging failed", 4)
  distribution = new_distribution

  tolerance = 50.0_dp * epsilon(1.0_dp) * &
    max(1.0_dp, config%final_time)
  do while (sparse_solution%time < config%final_time - tolerance)
    if (sparse_solution%steps >= config%maximum_steps) &
      call abort_run("Sparse MPI AMR step limit reached", 5)
    if (config%transport_enabled) then
      call sparse_patch_tree_reactive_timestep_1d( &
        species, config, distribution, sparse_solution, dt, ok, transport)
    else
      call sparse_patch_tree_reactive_timestep_1d( &
        species, config, distribution, sparse_solution, dt, ok)
    end if
    if (.not. ok .or. dt <= 0.0_dp) &
      call abort_run("Sparse MPI AMR timestep failed", 5)
    dt = min(dt, config%final_time - sparse_solution%time)
    if (config%transport_enabled) then
      call advance_sparse_patch_tree_reactive_1d( &
        species, reactions, config, dt, distribution, sparse_solution, &
        ok, transport)
    else
      call advance_sparse_patch_tree_reactive_1d( &
        species, reactions, config, dt, distribution, sparse_solution, ok)
    end if
    if (.not. ok) call abort_run("Sparse MPI AMR advance failed", 5)
    if (mod(sparse_solution%steps, config%amr_regrid_interval) == 0) then
      call regrid_tagged_sparse_patch_tree_reactive_1d( &
        species, config, distribution, sparse_solution, new_distribution, &
        changed, tagged_cells, transferred_cells, ok)
      if (.not. ok) call abort_run("Sparse MPI AMR regrid failed", 5)
      distribution = new_distribution
    end if
  end do
  sparse_solution%time = config%final_time

  call materialize_owned_patch_tree_reactive_1d( &
    distribution, sparse_solution, replicated_solution, ok)
  if (.not. ok) call abort_run("Final sparse AMR gather failed", 6)
  call patch_tree_reactive_integrals_1d(replicated_solution, final_all, ok)
  if (.not. ok) call abort_run("Final AMR integral failed", 6)
  final_integrals = final_all([irho, imx, imy, imz, iet])
  conservation_error = abs(final_integrals - initial_integrals) / &
    max(1.0_dp, abs(initial_integrals))

  output_ok = .true.
  if (rank == 0) then
    call write_patch_tree_reactive_1d_csv( &
      trim(output_path), species, replicated_solution, output_ok)
  end if
  call MPI_Bcast( &
    output_ok, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  if (ierr /= MPI_SUCCESS .or. .not. output_ok) &
    call abort_run("Sparse MPI AMR output failed", 6)

  if (rank == 0) then
    write(*, '(a)') "PeleF " // pelef_version // " sparse MPI AMR reactive 1D"
    write(*, '(a,i0)') "MPI ranks: ", nranks
    write(*, '(a,i0)') "Coarse cells: ", config%nx
    write(*, '(a,i0)') "Active AMR levels: ", &
      replicated_solution%level_count()
    write(*, '(a,i0)') "Active patches: ", &
      sum(distribution%rank_patch_counts)
    write(*, '(a,i0)') "Completed coarse steps: ", sparse_solution%steps
    write(*, '(a,i0)') "Regrid evaluations: ", &
      sparse_solution%regrid_evaluations
    write(*, '(a,i0)') "Hierarchy changes: ", sparse_solution%regrids
    write(*, '(a,i0)') "MPI work exponent: ", &
      config%amr_mpi_work_exponent
    write(*, '(a,es24.16)') "Final time: ", sparse_solution%time
    write(*, '(a,es24.16)') "Maximum conservation error: ", &
      maxval(conservation_error)
    write(*, '(a,1x,a)') "Output:", trim(output_path)
  end if

  call MPI_Finalize(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Finalize failed"

contains

  subroutine abort_run(reason, code)
    character(len=*), intent(in) :: reason
    integer, intent(in) :: code

    integer :: abort_ierr

    if (rank == 0) write(error_unit, '(a)') trim(reason)
    call MPI_Abort(MPI_COMM_WORLD, code, abort_ierr)
    error stop code
  end subroutine abort_run

end program pelef_mpi_amr_reactive_1d
