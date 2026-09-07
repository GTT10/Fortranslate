module mpi_amr_reactive_1d_application_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: error_unit
  use mpi_f08
  use precision_mod, only: dp
  use constants_mod, only: pelef_version
  use state_indices_mod, only: irho, imx, imy, imz, iet
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use simulation_config_reactive_1d_mod, only: reactive_1d_config
  use reactive_1d_mod, only: reactive_nvar
  use amr_patch_tree_1d_mod, only: amr_patch_level_plan_1d
  use amr_patch_tree_reactive_1d_mod, only: &
    amr_patch_tree_reactive_solution_1d, &
    initialize_patch_tree_reactive_1d, patch_tree_reactive_integrals_1d, &
    write_patch_tree_reactive_1d_checkpoint, &
    read_patch_tree_reactive_1d_checkpoint, &
    write_patch_tree_reactive_1d_selected_checkpoint, &
    read_patch_tree_reactive_1d_selected_checkpoint, &
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
  private

  public :: run_mpi_amr_reactive_1d_application

contains

  subroutine run_mpi_amr_reactive_1d_application( &
      comm, application_label, config, species, reactions, transport, &
      output_path, bundle_sha256, base_mole_fractions, chemistry_integrator)
    type(MPI_Comm), intent(in) :: comm
    character(len=*), intent(in) :: application_label, output_path
    type(reactive_1d_config), intent(in) :: config
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    character(len=*), intent(in), optional :: bundle_sha256
    real(dp), intent(in), optional :: base_mole_fractions(:)
    character(len=*), intent(in), optional :: chemistry_integrator

    type(amr_patch_level_plan_1d), allocatable :: empty_plans(:)
    type(amr_patch_tree_reactive_solution_1d) :: root_solution
    type(amr_patch_tree_reactive_solution_1d) :: replicated_solution
    type(amr_patch_tree_reactive_solution_1d) :: checkpoint_solution
    type(mpi_amr_sparse_reactive_solution_1d) :: sparse_solution
    type(mpi_amr_patch_distribution_1d) :: distribution, new_distribution
    real(dp), allocatable :: initial_all(:), final_all(:)
    real(dp), allocatable :: selected_composition(:)
    real(dp) :: initial_integrals(5), final_integrals(5)
    real(dp) :: conservation_error(5), dt, tolerance
    character(len=64) :: selected_bundle_sha256
    character(len=32) :: selected_chemistry_integrator
    logical :: ok, changed, output_ok, restart_run, selected_context
    logical :: stopped_after_checkpoint
    integer :: ierr, rank, nranks, tagged_cells, transferred_cells
    integer :: last_checkpoint_step, nvar

    call MPI_Comm_rank(comm, rank, ierr)
    if (ierr /= MPI_SUCCESS) error stop "MPI_Comm_rank failed"
    call MPI_Comm_size(comm, nranks, ierr)
    if (ierr /= MPI_SUCCESS) call abort_run("MPI_Comm_size failed", 2)

    selected_context = present(bundle_sha256) .or. &
      present(base_mole_fractions) .or. present(chemistry_integrator)
    if (selected_context .and. &
        (.not. present(bundle_sha256) .or. &
         .not. present(base_mole_fractions) .or. &
         .not. present(chemistry_integrator))) then
      call abort_run("Selected sparse MPI AMR context is incomplete", 2)
    end if
    if (selected_context) then
      selected_bundle_sha256 = bundle_sha256
      selected_chemistry_integrator = chemistry_integrator
      selected_composition = base_mole_fractions
    end if
    if (.not. config%amr_enabled .or. config%amr_multipatch_enabled) then
      call abort_run("Sparse MPI AMR patch-tree configuration is invalid", 2)
    end if
    nvar = reactive_nvar(size(species))
    if (nvar < 1 .or. size(reactions) < 1 .or. &
        size(transport) /= size(species) .or. len_trim(output_path) == 0) then
      call abort_run("Sparse MPI AMR application data are incompatible", 2)
    end if
    if (selected_context) then
      call selected_context_has_rank_consensus(ok)
      if (.not. ok) &
        call abort_consensus_run( &
          "Selected sparse MPI AMR rank context mismatch", 3)
    end if

    restart_run = len_trim(config%restart_file) > 0
    if (restart_run) then
      if (selected_context) then
        call read_patch_tree_reactive_1d_selected_checkpoint( &
          config%restart_file, species, config, selected_bundle_sha256, &
          selected_chemistry_integrator, selected_composition, &
          root_solution, initial_all, ok)
      else
        call read_patch_tree_reactive_1d_checkpoint( &
          config%restart_file, species, config, root_solution, ok)
      end if
      if (.not. ok) call abort_run("Sparse MPI AMR restart read failed", 4)
    else
      allocate(empty_plans(0))
      if (selected_context) then
        call initialize_patch_tree_reactive_1d( &
          species, config, empty_plans, root_solution, ok, &
          selected_composition)
      else
        call initialize_patch_tree_reactive_1d( &
          species, config, empty_plans, root_solution, ok)
      end if
      if (.not. ok) &
        call abort_run("Sparse MPI AMR root initialization failed", 4)
    end if
    tolerance = 50.0_dp * epsilon(1.0_dp) * &
      max(1.0_dp, config%final_time)
    if (root_solution%time > config%final_time + tolerance) &
      call abort_run("Restart time exceeds configured final_time", 4)
    if (.not. allocated(initial_all)) then
      allocate(initial_all(nvar))
      call patch_tree_reactive_integrals_1d(root_solution, initial_all, ok)
      if (.not. ok) call abort_run("Initial AMR integral failed", 4)
    end if
    allocate(final_all(nvar))
    initial_integrals = initial_all([irho, imx, imy, imz, iet])

    call initialize_mpi_amr_patch_distribution_1d( &
      root_solution%hierarchy, comm, distribution, ok, &
      config%amr_mpi_work_exponent)
    if (.not. ok) call abort_run("Initial sparse AMR ownership failed", 4)
    call scatter_owned_patch_tree_reactive_1d( &
      distribution, root_solution, sparse_solution, ok)
    if (.not. ok) call abort_run("Initial sparse AMR scatter failed", 4)
    if (.not. restart_run) then
      call regrid_tagged_sparse_patch_tree_reactive_1d( &
        species, config, distribution, sparse_solution, new_distribution, &
        changed, tagged_cells, transferred_cells, ok)
      if (.not. ok) call abort_run("Initial sparse AMR tagging failed", 4)
      distribution = new_distribution
    end if

    stopped_after_checkpoint = .false.
    last_checkpoint_step = -1
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
          ok, transport, chemistry_integrator=chemistry_integrator)
      else
        call advance_sparse_patch_tree_reactive_1d( &
          species, reactions, config, dt, distribution, sparse_solution, ok, &
          chemistry_integrator=chemistry_integrator)
      end if
      if (.not. ok) call abort_run("Sparse MPI AMR advance failed", 5)
      if (mod(sparse_solution%steps, config%amr_regrid_interval) == 0) then
        call regrid_tagged_sparse_patch_tree_reactive_1d( &
          species, config, distribution, sparse_solution, new_distribution, &
          changed, tagged_cells, transferred_cells, ok)
        if (.not. ok) call abort_run("Sparse MPI AMR regrid failed", 5)
        distribution = new_distribution
      end if
      if (config%checkpoint_interval > 0) then
        if (mod(sparse_solution%steps, config%checkpoint_interval) == 0) then
          call write_sparse_checkpoint(ok)
          if (.not. ok) call abort_run("Sparse MPI AMR checkpoint failed", 6)
          last_checkpoint_step = sparse_solution%steps
          if (config%checkpoint_stop_after_write) then
            stopped_after_checkpoint = .true.
            exit
          end if
        end if
      end if
    end do
    if (.not. stopped_after_checkpoint) &
      sparse_solution%time = config%final_time

    call materialize_owned_patch_tree_reactive_1d( &
      distribution, sparse_solution, replicated_solution, ok)
    if (.not. ok) call abort_run("Final sparse AMR gather failed", 6)
    if (len_trim(config%checkpoint_file) > 0 .and. &
        last_checkpoint_step /= sparse_solution%steps) then
      output_ok = .true.
      if (rank == 0) call write_materialized_checkpoint( &
        replicated_solution, output_ok)
      call MPI_Bcast(output_ok, 1, MPI_LOGICAL, 0, comm, ierr)
      if (ierr /= MPI_SUCCESS .or. .not. output_ok) &
        call abort_run("Final sparse MPI AMR checkpoint failed", 6)
    end if
    call patch_tree_reactive_integrals_1d( &
      replicated_solution, final_all, ok)
    if (.not. ok) call abort_run("Final AMR integral failed", 6)
    final_integrals = final_all([irho, imx, imy, imz, iet])
    conservation_error = abs(final_integrals - initial_integrals) / &
      max(1.0_dp, abs(initial_integrals))

    output_ok = .true.
    if (rank == 0) then
      call write_patch_tree_reactive_1d_csv( &
        trim(output_path), species, replicated_solution, output_ok)
    end if
    call MPI_Bcast(output_ok, 1, MPI_LOGICAL, 0, comm, ierr)
    if (ierr /= MPI_SUCCESS .or. .not. output_ok) &
      call abort_run("Sparse MPI AMR output failed", 6)

    if (rank == 0) then
      write(*, '(a)') "PeleF " // pelef_version // " " // &
        trim(application_label)
      if (selected_context) then
        write(*, '(a,1x,a)') &
          "Bundle SHA-256:", trim(selected_bundle_sha256)
        write(*, '(a,i0)') "Species: ", size(species)
        write(*, '(a,i0)') "Reactions: ", size(reactions)
        write(*, '(a,1x,a)') "Chemistry integrator:", &
          trim(selected_chemistry_integrator)
      end if
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
      write(*, '(a,l2)') "Restarted: ", restart_run
      write(*, '(a,l2)') "Stopped after checkpoint: ", &
        stopped_after_checkpoint
      write(*, '(a,es24.16)') "Final time: ", sparse_solution%time
      write(*, '(a,es24.16)') "Maximum conservation error: ", &
        maxval(conservation_error)
      write(*, '(a,1x,a)') "Output:", trim(output_path)
      if (len_trim(config%checkpoint_file) > 0) &
        write(*, '(a,1x,a)') "Checkpoint:", trim(config%checkpoint_file)
    end if

  contains

    subroutine selected_context_has_rank_consensus(context_ok)
      logical, intent(out) :: context_ok

      character(len=64) :: root_bundle
      character(len=32) :: root_integrator
      real(dp), allocatable :: root_composition(:)
      logical :: local_ok, global_ok
      integer :: context_ierr, index

      local_ok = len_trim(selected_bundle_sha256) == 64 .and. &
        size(selected_composition) == size(species) .and. &
        (trim(selected_chemistry_integrator) == "explicit" .or. &
         trim(selected_chemistry_integrator) == "implicit")
      if (local_ok) then
        do index = 1, 64
          select case (selected_bundle_sha256(index:index))
          case ('0':'9', 'a':'f', 'A':'F')
          case default
            local_ok = .false.
          end select
        end do
      end if
      if (local_ok) then
        local_ok = all(ieee_is_finite(selected_composition))
      end if
      if (local_ok) then
        local_ok = minval(selected_composition) >= 0.0_dp .and. &
          abs(sum(selected_composition) - 1.0_dp) <= 5.0e-10_dp
      end if
      call MPI_Allreduce( &
        local_ok, global_ok, 1, MPI_LOGICAL, MPI_LAND, comm, context_ierr)
      if (context_ierr /= MPI_SUCCESS .or. .not. global_ok) then
        context_ok = .false.
        return
      end if

      root_bundle = selected_bundle_sha256
      root_integrator = selected_chemistry_integrator
      allocate(root_composition(size(species)))
      root_composition = selected_composition
      call MPI_Bcast( &
        root_bundle, len(root_bundle), MPI_CHARACTER, 0, comm, context_ierr)
      if (context_ierr /= MPI_SUCCESS) then
        context_ok = .false.
        return
      end if
      call MPI_Bcast( &
        root_integrator, len(root_integrator), MPI_CHARACTER, 0, comm, &
        context_ierr)
      if (context_ierr /= MPI_SUCCESS) then
        context_ok = .false.
        return
      end if
      call MPI_Bcast( &
        root_composition, size(root_composition), MPI_DOUBLE_PRECISION, 0, &
        comm, context_ierr)
      if (context_ierr /= MPI_SUCCESS) then
        context_ok = .false.
        return
      end if
      local_ok = trim(root_bundle) == trim(selected_bundle_sha256) .and. &
        trim(root_integrator) == trim(selected_chemistry_integrator) .and. &
        .not. any(abs(root_composition - selected_composition) > 0.0_dp)
      call MPI_Allreduce( &
        local_ok, global_ok, 1, MPI_LOGICAL, MPI_LAND, comm, context_ierr)
      context_ok = context_ierr == MPI_SUCCESS .and. global_ok
    end subroutine selected_context_has_rank_consensus

    subroutine write_materialized_checkpoint(solution, checkpoint_ok)
      type(amr_patch_tree_reactive_solution_1d), intent(in) :: solution
      logical, intent(out) :: checkpoint_ok

      if (selected_context) then
        call write_patch_tree_reactive_1d_selected_checkpoint( &
          config%checkpoint_file, species, solution, &
          selected_bundle_sha256, selected_chemistry_integrator, &
          selected_composition, initial_all, checkpoint_ok)
      else
        call write_patch_tree_reactive_1d_checkpoint( &
          config%checkpoint_file, species, solution, checkpoint_ok)
      end if
    end subroutine write_materialized_checkpoint

    subroutine write_sparse_checkpoint(checkpoint_ok)
      logical, intent(out) :: checkpoint_ok

      call materialize_owned_patch_tree_reactive_1d( &
        distribution, sparse_solution, checkpoint_solution, checkpoint_ok)
      if (.not. checkpoint_ok) return
      if (rank == 0) call write_materialized_checkpoint( &
        checkpoint_solution, checkpoint_ok)
      call MPI_Bcast(checkpoint_ok, 1, MPI_LOGICAL, 0, comm, ierr)
      checkpoint_ok = ierr == MPI_SUCCESS .and. checkpoint_ok
    end subroutine write_sparse_checkpoint

    subroutine abort_consensus_run(reason, code)
      character(len=*), intent(in) :: reason
      integer, intent(in) :: code

      integer :: abort_ierr, sync_ierr

      if (rank == 0) then
        write(error_unit, '(a)') trim(reason)
        flush(error_unit)
      end if
      call MPI_Barrier(comm, sync_ierr)
      call MPI_Abort(comm, code, abort_ierr)
      error stop code
    end subroutine abort_consensus_run

    subroutine abort_run(reason, code)
      character(len=*), intent(in) :: reason
      integer, intent(in) :: code

      integer :: abort_ierr

      if (rank == 0) then
        write(error_unit, '(a)') trim(reason)
        flush(error_unit)
      end if
      call MPI_Abort(comm, code, abort_ierr)
      error stop code
    end subroutine abort_run

  end subroutine run_mpi_amr_reactive_1d_application

end module mpi_amr_reactive_1d_application_mod
