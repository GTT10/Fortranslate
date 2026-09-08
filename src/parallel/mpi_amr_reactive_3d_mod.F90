module mpi_amr_reactive_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: int64
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_conserved_to_primitive, &
    reactive_riemann_flux_x
  use reactive_2d_mod, only: reactive_riemann_flux_y
  use reactive_directional_flux_3d_mod, only: reactive_riemann_flux_z
  use reactive_3d_mod, only: &
    recover_reactive_temperatures_3d, &
    compute_reactive_plm_slab_face_fluxes_3d
  use amr_hierarchy_3d_mod, only: amr_patch_3d, average_down_3d
  use amr_reactive_3d_mod, only: &
    accumulate_fine_interface_fluxes_3d, reflux_coarse_3d, &
    compute_fine_patch_plm_slab_face_fluxes_3d
  implicit none
  private

  integer, parameter :: maximum_species_count = 32
  integer, parameter :: nasa7_real_count = 18

  type, public :: mpi_amr_slab_distribution_3d
    type(MPI_Comm) :: comm = MPI_COMM_NULL
    integer :: rank = -1
    integer :: nranks = 0
    integer :: coarse_first = 0
    integer :: coarse_last = -1
    integer :: fine_first = 0
    integer :: fine_last = -1
    integer, allocatable :: coarse_counts(:)
    integer, allocatable :: coarse_displacements(:)
    integer, allocatable :: fine_counts(:)
    integer, allocatable :: fine_displacements(:)
  contains
    procedure :: is_valid => mpi_amr_slab_distribution_is_valid_3d
  end type mpi_amr_slab_distribution_3d

  public :: initialize_mpi_amr_slab_distribution_3d
  public :: compute_mpi_amr_reactive_cfl_timestep_3d
  public :: advance_mpi_amr_reactive_hydro_3d
  public :: broadcast_mpi_amr_reactive_hierarchy_3d

contains

  pure logical function mpi_amr_slab_distribution_is_valid_3d( &
      self, patch) result(valid)
    class(mpi_amr_slab_distribution_3d), intent(in) :: self
    type(amr_patch_3d), intent(in) :: patch

    valid = patch%is_strictly_interior() .and. self%nranks >= 1 .and. &
      self%rank >= 0 .and. self%rank < self%nranks .and. &
      allocated(self%coarse_counts) .and. &
      allocated(self%coarse_displacements) .and. &
      allocated(self%fine_counts) .and. &
      allocated(self%fine_displacements)
    if (.not. valid) return
    valid = size(self%coarse_counts) == self%nranks .and. &
      size(self%coarse_displacements) == self%nranks .and. &
      size(self%fine_counts) == self%nranks .and. &
      size(self%fine_displacements) == self%nranks .and. &
      all(self%coarse_counts >= 1) .and. all(self%fine_counts >= 1) .and. &
      sum(self%coarse_counts) == patch%coarse_nx .and. &
      sum(self%fine_counts) == patch%fine_nx()
    if (.not. valid) return
    valid = self%coarse_displacements(1) == 0 .and. &
      self%fine_displacements(1) == 0
    if (.not. valid) return
    if (self%nranks > 1) then
      valid = all(self%coarse_displacements(2:) == &
        self%coarse_displacements(:self%nranks - 1) + &
          self%coarse_counts(:self%nranks - 1)) .and. &
        all(self%fine_displacements(2:) == &
          self%fine_displacements(:self%nranks - 1) + &
            self%fine_counts(:self%nranks - 1))
      if (.not. valid) return
    end if
    valid = self%coarse_first == &
      self%coarse_displacements(self%rank + 1) + 1 .and. &
      self%coarse_last == self%coarse_first + &
        self%coarse_counts(self%rank + 1) - 1 .and. &
      self%fine_first == self%fine_displacements(self%rank + 1) + 1 .and. &
      self%fine_last == self%fine_first + &
        self%fine_counts(self%rank + 1) - 1
  end function mpi_amr_slab_distribution_is_valid_3d

  subroutine initialize_mpi_amr_slab_distribution_3d( &
      patch, comm, distribution, ok)
    type(amr_patch_3d), intent(in) :: patch
    type(MPI_Comm), intent(in) :: comm
    type(mpi_amr_slab_distribution_3d), intent(out) :: distribution
    logical, intent(out) :: ok

    integer :: ierr

    ok = .false.
    distribution%comm = comm
    call MPI_Comm_rank(comm, distribution%rank, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Comm_size(comm, distribution%nranks, ierr)
    if (ierr /= MPI_SUCCESS) return
    if (distribution%nranks < 1) return
    if (.not. patch%is_strictly_interior() .or. &
        distribution%nranks > patch%coarse_nx .or. &
        distribution%nranks > patch%fine_nx()) return
    allocate(distribution%coarse_counts(distribution%nranks))
    allocate(distribution%coarse_displacements(distribution%nranks))
    allocate(distribution%fine_counts(distribution%nranks))
    allocate(distribution%fine_displacements(distribution%nranks))
    call build_slab_counts( &
      patch%coarse_nx, distribution%nranks, &
      distribution%coarse_counts, distribution%coarse_displacements)
    call build_slab_counts( &
      patch%fine_nx(), distribution%nranks, &
      distribution%fine_counts, distribution%fine_displacements)
    distribution%coarse_first = &
      distribution%coarse_displacements(distribution%rank + 1) + 1
    distribution%coarse_last = distribution%coarse_first + &
      distribution%coarse_counts(distribution%rank + 1) - 1
    distribution%fine_first = &
      distribution%fine_displacements(distribution%rank + 1) + 1
    distribution%fine_last = distribution%fine_first + &
      distribution%fine_counts(distribution%rank + 1) - 1
    ok = distribution%is_valid(patch)
  end subroutine initialize_mpi_amr_slab_distribution_3d

  pure subroutine build_slab_counts(global_count, nranks, counts, displacements)
    integer, intent(in) :: global_count, nranks
    integer, intent(out) :: counts(:), displacements(:)

    integer :: rank_index, base_count, remainder

    base_count = global_count / nranks
    remainder = modulo(global_count, nranks)
    displacements(1) = 0
    do rank_index = 1, nranks
      counts(rank_index) = base_count
      if (rank_index <= remainder) counts(rank_index) = &
        counts(rank_index) + 1
    end do
    do rank_index = 2, nranks
      displacements(rank_index) = displacements(rank_index - 1) + &
        counts(rank_index - 1)
    end do
  end subroutine build_slab_counts

  subroutine compute_mpi_amr_reactive_cfl_timestep_3d( &
      distribution, species, patch, coarse_state, coarse_temperature, &
      fine_state, fine_temperature, dx, dy, dz, cfl, dt, ok)
    type(mpi_amr_slab_distribution_3d), intent(in) :: distribution
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: coarse_state(:, :, :, :)
    real(dp), intent(in) :: coarse_temperature(:, :, :)
    real(dp), intent(in) :: fine_state(:, :, :, :)
    real(dp), intent(in) :: fine_temperature(:, :, :)
    real(dp), intent(in) :: dx, dy, dz, cfl
    real(dp), intent(out) :: dt
    logical, intent(out) :: ok

    real(dp) :: coarse_rate, fine_rate, ratio_real
    real(dp) :: metadata(4)
    logical :: local_ok

    dt = 0.0_dp
    metadata = [dx, dy, dz, cfl]
    local_ok = distribution%is_valid(patch) .and. &
      valid_mpi_amr_reactive_shapes( &
        species, patch, coarse_state, coarse_temperature, &
        fine_state, fine_temperature)
    if (local_ok) local_ok = all(ieee_is_finite(metadata))
    if (local_ok) local_ok = &
      dx > 0.0_dp .and. dy > 0.0_dp .and. dz > 0.0_dp .and. &
      cfl > 0.0_dp .and. cfl <= 1.0_dp
    call collective_metadata_matches( &
      distribution, species, patch, metadata, "", local_ok, ok)
    if (.not. ok) return
    call distributed_maximum_rate_3d( &
      distribution, species, coarse_state, coarse_temperature, &
      distribution%coarse_first, distribution%coarse_last, &
      dx, dy, dz, coarse_rate, ok)
    if (.not. ok) return
    ratio_real = real(patch%refinement_ratio, dp)
    call distributed_maximum_rate_3d( &
      distribution, species, fine_state, fine_temperature, &
      distribution%fine_first, distribution%fine_last, &
      dx / ratio_real, dy / ratio_real, dz / ratio_real, fine_rate, ok)
    if (.not. ok) return
    dt = min(cfl / coarse_rate, ratio_real * cfl / fine_rate)
    ok = ieee_is_finite(dt)
    if (ok) ok = dt > 0.0_dp
    call collective_logical_and(distribution%comm, ok, local_ok)
    ok = local_ok
    if (.not. ok) dt = 0.0_dp
  end subroutine compute_mpi_amr_reactive_cfl_timestep_3d

  subroutine distributed_maximum_rate_3d( &
      distribution, species, state, temperature, first_i, last_i, &
      dx, dy, dz, maximum_rate, ok)
    type(mpi_amr_slab_distribution_3d), intent(in) :: distribution
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: first_i, last_i
    real(dp), intent(in) :: dx, dy, dz
    real(dp), intent(out) :: maximum_rate
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:)
    real(dp) :: local_maximum, cell_temperature, sound_speed, rate
    logical :: cell_ok, local_ok, global_ok
    integer :: i, j, k, ierr

    allocate(primitive(reactive_nprim(size(species))))
    local_maximum = 0.0_dp
    local_ok = first_i >= 1 .and. last_i <= size(state, 2) .and. &
      last_i >= first_i
    if (local_ok) then
      outer: do k = 1, size(state, 4)
        do j = 1, size(state, 3)
          do i = first_i, last_i
            call reactive_conserved_to_primitive( &
              species, state(:, i, j, k), temperature(i, j, k), &
              primitive, cell_temperature, sound_speed, cell_ok)
            if (.not. cell_ok) then
              local_ok = .false.
              exit outer
            end if
            rate = (abs(primitive(2)) + sound_speed) / dx + &
              (abs(primitive(3)) + sound_speed) / dy + &
              (abs(primitive(4)) + sound_speed) / dz
            local_maximum = max(local_maximum, rate)
          end do
        end do
      end do outer
    end if
    call collective_logical_and(distribution%comm, local_ok, global_ok)
    if (.not. global_ok) then
      maximum_rate = 0.0_dp
      ok = .false.
      return
    end if
    call MPI_Allreduce( &
      local_maximum, maximum_rate, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
      distribution%comm, ierr)
    ok = ierr == MPI_SUCCESS
    if (ok) ok = ieee_is_finite(maximum_rate)
    if (ok) ok = maximum_rate > 0.0_dp
  end subroutine distributed_maximum_rate_3d

  subroutine advance_mpi_amr_reactive_hydro_3d( &
      distribution, species, patch, coarse_state, coarse_temperature, &
      fine_state, fine_temperature, dx, dy, dz, dt, riemann_solver, &
      maximum_reflux_correction, ok, reconstruction, limiter)
    type(mpi_amr_slab_distribution_3d), intent(in) :: distribution
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(inout) :: coarse_state(:, :, :, :)
    real(dp), intent(inout) :: coarse_temperature(:, :, :)
    real(dp), intent(inout) :: fine_state(:, :, :, :)
    real(dp), intent(inout) :: fine_temperature(:, :, :)
    real(dp), intent(in) :: dx, dy, dz, dt
    character(len=*), intent(in) :: riemann_solver
    real(dp), intent(out) :: maximum_reflux_correction
    logical, intent(out) :: ok
    character(len=*), intent(in), optional :: reconstruction, limiter

    real(dp), allocatable :: coarse_start(:, :, :, :)
    real(dp), allocatable :: coarse_candidate(:, :, :, :)
    real(dp), allocatable :: fine_candidate(:, :, :, :)
    real(dp), allocatable :: coarse_start_temperature(:, :, :)
    real(dp), allocatable :: coarse_candidate_temperature(:, :, :)
    real(dp), allocatable :: fine_candidate_temperature(:, :, :)
    real(dp), allocatable :: synchronized_temperature(:, :, :)
    real(dp), allocatable :: coarse_flux_x(:, :, :, :)
    real(dp), allocatable :: coarse_flux_y(:, :, :, :)
    real(dp), allocatable :: coarse_flux_z(:, :, :, :)
    real(dp), allocatable :: fine_flux_x(:, :, :, :)
    real(dp), allocatable :: fine_flux_y(:, :, :, :)
    real(dp), allocatable :: fine_flux_z(:, :, :, :)
    real(dp), allocatable :: fine_x_lower(:, :, :), fine_x_upper(:, :, :)
    real(dp), allocatable :: fine_y_lower(:, :, :), fine_y_upper(:, :, :)
    real(dp), allocatable :: fine_z_lower(:, :, :), fine_z_upper(:, :, :)
    real(dp) :: fine_dt, alpha_start, alpha_end, metadata(4)
    real(dp) :: local_reflux, global_reflux
    logical :: local_ok, global_ok
    integer :: nvar, ratio, substep, ierr
    integer :: covered_nx, covered_ny, covered_nz
    character(len=32) :: selected_reconstruction, selected_limiter
    character(len=256) :: algorithm_contract

    maximum_reflux_correction = 0.0_dp
    ok = .false.
    selected_reconstruction = "pcm"
    selected_limiter = "mc"
    if (present(reconstruction)) selected_reconstruction = &
      trim(reconstruction)
    if (present(limiter)) selected_limiter = trim(limiter)
    algorithm_contract = trim(riemann_solver) // "|" // &
      trim(selected_reconstruction) // "|" // trim(selected_limiter)
    metadata = [dx, dy, dz, dt]
    local_ok = distribution%is_valid(patch) .and. &
      valid_mpi_amr_reactive_shapes( &
        species, patch, coarse_state, coarse_temperature, &
        fine_state, fine_temperature)
    if (local_ok) local_ok = all(ieee_is_finite(metadata))
    if (local_ok) local_ok = &
      dx > 0.0_dp .and. dy > 0.0_dp .and. dz > 0.0_dp .and. &
      dt > 0.0_dp
    if (local_ok) local_ok = &
      trim(selected_reconstruction) == "pcm" .or. &
      trim(selected_reconstruction) == "characteristic_plm"
    if (local_ok) local_ok = &
      trim(selected_limiter) == "minmod" .or. &
      trim(selected_limiter) == "mc"
    call collective_metadata_matches( &
      distribution, species, patch, metadata, algorithm_contract, local_ok, ok)
    if (.not. ok) return

    nvar = reactive_nvar(size(species))
    ratio = patch%refinement_ratio
    covered_nx = patch%coarse_i_upper - patch%coarse_i_lower + 1
    covered_ny = patch%coarse_j_upper - patch%coarse_j_lower + 1
    covered_nz = patch%coarse_k_upper - patch%coarse_k_lower + 1
    coarse_start = coarse_state
    coarse_candidate = coarse_state
    fine_candidate = fine_state
    allocate(coarse_start_temperature, mold=coarse_temperature)
    coarse_candidate_temperature = coarse_temperature
    fine_candidate_temperature = fine_temperature
    call recover_reactive_temperatures_3d( &
      species, coarse_start, coarse_temperature, patch%coarse_nx, &
      patch%coarse_ny, patch%coarse_nz, coarse_start_temperature, local_ok)
    call collective_logical_and(distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return

    allocate(coarse_flux_x, mold=coarse_state)
    allocate(coarse_flux_y, mold=coarse_state)
    allocate(coarse_flux_z, mold=coarse_state)
    call distributed_periodic_ssprk2_3d( &
      distribution, species, distribution%coarse_first, &
      distribution%coarse_last, coarse_candidate, &
      coarse_candidate_temperature, dx, dy, dz, dt, riemann_solver, &
      selected_reconstruction, selected_limiter, coarse_flux_x, &
      coarse_flux_y, coarse_flux_z, ok)
    if (.not. ok) return

    allocate(fine_flux_x(nvar, 0:patch%fine_nx(), &
      patch%fine_ny(), patch%fine_nz()))
    allocate(fine_flux_y(nvar, patch%fine_nx(), &
      0:patch%fine_ny(), patch%fine_nz()))
    allocate(fine_flux_z(nvar, patch%fine_nx(), &
      patch%fine_ny(), 0:patch%fine_nz()))
    allocate(fine_x_lower(nvar, covered_ny, covered_nz), &
      fine_x_upper(nvar, covered_ny, covered_nz))
    allocate(fine_y_lower(nvar, covered_nx, covered_nz), &
      fine_y_upper(nvar, covered_nx, covered_nz))
    allocate(fine_z_lower(nvar, covered_nx, covered_ny), &
      fine_z_upper(nvar, covered_nx, covered_ny))
    fine_x_lower = 0.0_dp
    fine_x_upper = 0.0_dp
    fine_y_lower = 0.0_dp
    fine_y_upper = 0.0_dp
    fine_z_lower = 0.0_dp
    fine_z_upper = 0.0_dp

    fine_dt = dt / real(ratio, dp)
    do substep = 1, ratio
      alpha_start = real(substep - 1, dp) / real(ratio, dp)
      alpha_end = real(substep, dp) / real(ratio, dp)
      call distributed_fine_patch_ssprk2_3d( &
        distribution, species, patch, coarse_start, &
        coarse_start_temperature, coarse_candidate, &
        coarse_candidate_temperature, fine_candidate, &
        fine_candidate_temperature, dx / real(ratio, dp), &
        dy / real(ratio, dp), dz / real(ratio, dp), fine_dt, &
        alpha_start, alpha_end, riemann_solver, selected_reconstruction, &
        selected_limiter, fine_flux_x, fine_flux_y, fine_flux_z, ok)
      if (.not. ok) return
      call accumulate_fine_interface_fluxes_3d( &
        patch, fine_flux_x, fine_flux_y, fine_flux_z, &
        fine_x_lower, fine_x_upper, fine_y_lower, fine_y_upper, &
        fine_z_lower, fine_z_upper, local_ok)
      call collective_logical_and(distribution%comm, local_ok, global_ok)
      if (.not. global_ok) return
    end do
    fine_x_lower = fine_x_lower / real(ratio, dp)
    fine_x_upper = fine_x_upper / real(ratio, dp)
    fine_y_lower = fine_y_lower / real(ratio, dp)
    fine_y_upper = fine_y_upper / real(ratio, dp)
    fine_z_lower = fine_z_lower / real(ratio, dp)
    fine_z_upper = fine_z_upper / real(ratio, dp)

    call reflux_coarse_3d( &
      patch, coarse_candidate, coarse_flux_x, coarse_flux_y, coarse_flux_z, &
      fine_x_lower, fine_x_upper, fine_y_lower, fine_y_upper, &
      fine_z_lower, fine_z_upper, dx, dy, dz, dt, local_reflux, local_ok)
    call collective_logical_and(distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    call MPI_Allreduce( &
      local_reflux, global_reflux, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call average_down_3d(coarse_candidate, fine_candidate, patch, local_ok)
    call collective_logical_and(distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    allocate(synchronized_temperature, mold=coarse_temperature)
    call recover_reactive_temperatures_3d( &
      species, coarse_candidate, coarse_candidate_temperature, &
      patch%coarse_nx, patch%coarse_ny, patch%coarse_nz, &
      synchronized_temperature, local_ok)
    call collective_logical_and(distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return

    coarse_state = coarse_candidate
    coarse_temperature = synchronized_temperature
    fine_state = fine_candidate
    fine_temperature = fine_candidate_temperature
    maximum_reflux_correction = global_reflux
    ok = .true.
  end subroutine advance_mpi_amr_reactive_hydro_3d

  subroutine distributed_periodic_ssprk2_3d( &
      distribution, species, first_i, last_i, state, temperature, &
      dx, dy, dz, dt, riemann_solver, reconstruction, limiter, &
      face_flux_x, face_flux_y, face_flux_z, ok)
    type(mpi_amr_slab_distribution_3d), intent(in) :: distribution
    type(nasa7_species), intent(in) :: species(:)
    integer, intent(in) :: first_i, last_i
    real(dp), intent(inout) :: state(:, :, :, :), temperature(:, :, :)
    real(dp), intent(in) :: dx, dy, dz, dt
    character(len=*), intent(in) :: riemann_solver, reconstruction, limiter
    real(dp), intent(out) :: face_flux_x(:, :, :, :)
    real(dp), intent(out) :: face_flux_y(:, :, :, :)
    real(dp), intent(out) :: face_flux_z(:, :, :, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: old_state(:, :, :, :)
    real(dp), allocatable :: stage_state(:, :, :, :)
    real(dp), allocatable :: candidate_state(:, :, :, :)
    real(dp), allocatable :: old_temperature(:, :, :)
    real(dp), allocatable :: stage_temperature(:, :, :)
    real(dp), allocatable :: candidate_temperature(:, :, :)
    real(dp), allocatable :: first_flux_x(:, :, :, :)
    real(dp), allocatable :: first_flux_y(:, :, :, :)
    real(dp), allocatable :: first_flux_z(:, :, :, :)
    real(dp), allocatable :: second_flux_x(:, :, :, :)
    real(dp), allocatable :: second_flux_y(:, :, :, :)
    real(dp), allocatable :: second_flux_z(:, :, :, :)
    logical :: local_ok, global_ok

    face_flux_x = 0.0_dp
    face_flux_y = 0.0_dp
    face_flux_z = 0.0_dp
    ok = .false.
    old_state = state
    allocate(old_temperature, mold=temperature)
    call recover_reactive_temperatures_3d( &
      species, old_state, temperature, size(state, 2), size(state, 3), &
      size(state, 4), old_temperature, local_ok)
    call collective_logical_and(distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    allocate(first_flux_x, mold=state)
    allocate(first_flux_y, mold=state)
    allocate(first_flux_z, mold=state)
    call distributed_periodic_face_fluxes_3d( &
      distribution, species, first_i, last_i, old_state, old_temperature, &
      riemann_solver, reconstruction, limiter, first_flux_x, first_flux_y, &
      first_flux_z, ok)
    if (.not. ok) return
    call distributed_periodic_stage_3d( &
      distribution, first_i, last_i, old_state, old_state, &
      first_flux_x, first_flux_y, first_flux_z, dx, dy, dz, dt, &
      .false., stage_state, ok)
    if (.not. ok) return
    allocate(stage_temperature, mold=temperature)
    call recover_reactive_temperatures_3d( &
      species, stage_state, old_temperature, size(state, 2), &
      size(state, 3), size(state, 4), stage_temperature, local_ok)
    call collective_logical_and(distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return

    allocate(second_flux_x, mold=state)
    allocate(second_flux_y, mold=state)
    allocate(second_flux_z, mold=state)
    call distributed_periodic_face_fluxes_3d( &
      distribution, species, first_i, last_i, stage_state, &
      stage_temperature, riemann_solver, reconstruction, limiter, &
      second_flux_x, second_flux_y, second_flux_z, ok)
    if (.not. ok) return
    call distributed_periodic_stage_3d( &
      distribution, first_i, last_i, old_state, stage_state, &
      second_flux_x, second_flux_y, second_flux_z, dx, dy, dz, dt, &
      .true., candidate_state, ok)
    if (.not. ok) return
    allocate(candidate_temperature, mold=temperature)
    call recover_reactive_temperatures_3d( &
      species, candidate_state, stage_temperature, size(state, 2), &
      size(state, 3), size(state, 4), candidate_temperature, local_ok)
    call collective_logical_and(distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return

    face_flux_x = 0.5_dp * (first_flux_x + second_flux_x)
    face_flux_y = 0.5_dp * (first_flux_y + second_flux_y)
    face_flux_z = 0.5_dp * (first_flux_z + second_flux_z)
    state = candidate_state
    temperature = candidate_temperature
    ok = .true.
  end subroutine distributed_periodic_ssprk2_3d

  subroutine distributed_periodic_face_fluxes_3d( &
      distribution, species, first_i, last_i, state, temperature, &
      riemann_solver, reconstruction, limiter, flux_x, flux_y, flux_z, ok)
    type(mpi_amr_slab_distribution_3d), intent(in) :: distribution
    type(nasa7_species), intent(in) :: species(:)
    integer, intent(in) :: first_i, last_i
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    character(len=*), intent(in) :: riemann_solver, reconstruction, limiter
    real(dp), intent(out) :: flux_x(:, :, :, :)
    real(dp), intent(out) :: flux_y(:, :, :, :)
    real(dp), intent(out) :: flux_z(:, :, :, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: local_x(:, :, :, :)
    real(dp), allocatable :: local_y(:, :, :, :)
    real(dp), allocatable :: local_z(:, :, :, :)
    logical :: face_ok, local_ok, global_ok
    integer :: i, j, k, next_i, next_j, next_k

    flux_x = 0.0_dp
    flux_y = 0.0_dp
    flux_z = 0.0_dp
    allocate(local_x, mold=state)
    allocate(local_y, mold=state)
    allocate(local_z, mold=state)
    local_x = 0.0_dp
    local_y = 0.0_dp
    local_z = 0.0_dp
    local_ok = first_i >= 1 .and. last_i <= size(state, 2) .and. &
      last_i >= first_i .and. all(shape(flux_x) == shape(state)) .and. &
      all(shape(flux_y) == shape(state)) .and. &
      all(shape(flux_z) == shape(state)) .and. &
      (trim(reconstruction) == "pcm" .or. &
        trim(reconstruction) == "characteristic_plm") .and. &
      (trim(limiter) == "minmod" .or. trim(limiter) == "mc")
    if (local_ok) then
      select case (trim(reconstruction))
      case ("characteristic_plm")
        call compute_reactive_plm_slab_face_fluxes_3d( &
          species, state, temperature, size(state, 2), size(state, 3), &
          size(state, 4), first_i, last_i, limiter, riemann_solver, &
          local_x, local_y, local_z, face_ok)
        local_ok = face_ok
      case ("pcm")
        outer: do k = 1, size(state, 4)
          next_k = modulo(k, size(state, 4)) + 1
          do j = 1, size(state, 3)
            next_j = modulo(j, size(state, 3)) + 1
            do i = first_i, last_i
              next_i = modulo(i, size(state, 2)) + 1
              call reactive_riemann_flux_x( &
                species, state(:, i, j, k), state(:, next_i, j, k), &
                temperature(i, j, k), temperature(next_i, j, k), &
                riemann_solver, local_x(:, i, j, k), face_ok)
              if (.not. face_ok) then
                local_ok = .false.
                exit outer
              end if
              call reactive_riemann_flux_y( &
                species, state(:, i, j, k), state(:, i, next_j, k), &
                temperature(i, j, k), temperature(i, next_j, k), &
                riemann_solver, local_y(:, i, j, k), face_ok)
              if (.not. face_ok) then
                local_ok = .false.
                exit outer
              end if
              call reactive_riemann_flux_z( &
                species, state(:, i, j, k), state(:, i, j, next_k), &
                temperature(i, j, k), temperature(i, j, next_k), &
                riemann_solver, local_z(:, i, j, k), face_ok)
              if (.not. face_ok) then
                local_ok = .false.
                exit outer
              end if
            end do
          end do
        end do outer
      end select
    end if
    call collective_logical_and(distribution%comm, local_ok, global_ok)
    if (.not. global_ok) then
      ok = .false.
      return
    end if
    call collective_sum_rank4(distribution%comm, local_x, flux_x, ok)
    if (.not. ok) return
    call collective_sum_rank4(distribution%comm, local_y, flux_y, ok)
    if (.not. ok) return
    call collective_sum_rank4(distribution%comm, local_z, flux_z, ok)
    if (.not. ok) return
    local_ok = all(ieee_is_finite(flux_x)) .and. &
      all(ieee_is_finite(flux_y)) .and. all(ieee_is_finite(flux_z))
    call collective_logical_and(distribution%comm, local_ok, ok)
  end subroutine distributed_periodic_face_fluxes_3d

  subroutine distributed_periodic_stage_3d( &
      distribution, first_i, last_i, old_state, stage_base, &
      flux_x, flux_y, flux_z, dx, dy, dz, dt, second_stage, &
      candidate, ok)
    type(mpi_amr_slab_distribution_3d), intent(in) :: distribution
    integer, intent(in) :: first_i, last_i
    real(dp), intent(in) :: old_state(:, :, :, :)
    real(dp), intent(in) :: stage_base(:, :, :, :)
    real(dp), intent(in) :: flux_x(:, :, :, :)
    real(dp), intent(in) :: flux_y(:, :, :, :)
    real(dp), intent(in) :: flux_z(:, :, :, :)
    real(dp), intent(in) :: dx, dy, dz, dt
    logical, intent(in) :: second_stage
    real(dp), allocatable, intent(out) :: candidate(:, :, :, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: local_candidate(:, :, :, :), rhs(:)
    logical :: local_ok, global_ok
    integer :: i, j, k, previous_i, previous_j, previous_k

    allocate(local_candidate, mold=old_state)
    allocate(candidate, mold=old_state)
    allocate(rhs(size(old_state, 1)))
    local_candidate = 0.0_dp
    candidate = 0.0_dp
    local_ok = first_i >= 1 .and. last_i <= size(old_state, 2) .and. &
      last_i >= first_i .and. all(shape(stage_base) == shape(old_state))
    if (local_ok) then
      do k = 1, size(old_state, 4)
        previous_k = modulo(k - 2, size(old_state, 4)) + 1
        do j = 1, size(old_state, 3)
          previous_j = modulo(j - 2, size(old_state, 3)) + 1
          do i = first_i, last_i
            previous_i = modulo(i - 2, size(old_state, 2)) + 1
            rhs = -(flux_x(:, i, j, k) - &
                flux_x(:, previous_i, j, k)) / dx &
              -(flux_y(:, i, j, k) - &
                flux_y(:, i, previous_j, k)) / dy &
              -(flux_z(:, i, j, k) - &
                flux_z(:, i, j, previous_k)) / dz
            if (second_stage) then
              local_candidate(:, i, j, k) = 0.5_dp * &
                old_state(:, i, j, k) + 0.5_dp * &
                  (stage_base(:, i, j, k) + dt * rhs)
            else
              local_candidate(:, i, j, k) = &
                stage_base(:, i, j, k) + dt * rhs
            end if
          end do
        end do
      end do
      local_ok = all(ieee_is_finite( &
        local_candidate(:, first_i:last_i, :, :)))
    end if
    call collective_logical_and(distribution%comm, local_ok, global_ok)
    if (.not. global_ok) then
      ok = .false.
      return
    end if
    call collective_sum_rank4( &
      distribution%comm, local_candidate, candidate, ok)
  end subroutine distributed_periodic_stage_3d

  subroutine distributed_fine_patch_ssprk2_3d( &
      distribution, species, patch, coarse_start, &
      coarse_start_temperature, coarse_end, coarse_end_temperature, &
      state, temperature, dx, dy, dz, dt, alpha_start, alpha_end, &
      riemann_solver, reconstruction, limiter, face_flux_x, face_flux_y, &
      face_flux_z, ok)
    type(mpi_amr_slab_distribution_3d), intent(in) :: distribution
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: coarse_start(:, :, :, :)
    real(dp), intent(in) :: coarse_start_temperature(:, :, :)
    real(dp), intent(in) :: coarse_end(:, :, :, :)
    real(dp), intent(in) :: coarse_end_temperature(:, :, :)
    real(dp), intent(inout) :: state(:, :, :, :), temperature(:, :, :)
    real(dp), intent(in) :: dx, dy, dz, dt, alpha_start, alpha_end
    character(len=*), intent(in) :: riemann_solver, reconstruction, limiter
    real(dp), intent(out) :: face_flux_x(:, 0:, :, :)
    real(dp), intent(out) :: face_flux_y(:, :, 0:, :)
    real(dp), intent(out) :: face_flux_z(:, :, :, 0:)
    logical, intent(out) :: ok

    real(dp), allocatable :: old_state(:, :, :, :)
    real(dp), allocatable :: stage_state(:, :, :, :)
    real(dp), allocatable :: candidate_state(:, :, :, :)
    real(dp), allocatable :: old_temperature(:, :, :)
    real(dp), allocatable :: stage_temperature(:, :, :)
    real(dp), allocatable :: candidate_temperature(:, :, :)
    real(dp), allocatable :: first_flux_x(:, :, :, :)
    real(dp), allocatable :: first_flux_y(:, :, :, :)
    real(dp), allocatable :: first_flux_z(:, :, :, :)
    real(dp), allocatable :: second_flux_x(:, :, :, :)
    real(dp), allocatable :: second_flux_y(:, :, :, :)
    real(dp), allocatable :: second_flux_z(:, :, :, :)
    logical :: local_ok, global_ok
    integer :: nvar, nx, ny, nz

    face_flux_x = 0.0_dp
    face_flux_y = 0.0_dp
    face_flux_z = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    nx = patch%fine_nx()
    ny = patch%fine_ny()
    nz = patch%fine_nz()
    old_state = state
    allocate(old_temperature(nx, ny, nz))
    call recover_reactive_temperatures_3d( &
      species, old_state, temperature, nx, ny, nz, &
      old_temperature, local_ok)
    call collective_logical_and(distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    allocate(first_flux_x(nvar, 0:nx, ny, nz))
    allocate(first_flux_y(nvar, nx, 0:ny, nz))
    allocate(first_flux_z(nvar, nx, ny, 0:nz))
    call distributed_fine_face_fluxes_3d( &
      distribution, species, patch, coarse_start, &
      coarse_start_temperature, coarse_end, coarse_end_temperature, &
      old_state, old_temperature, alpha_start, riemann_solver, &
      reconstruction, limiter, first_flux_x, first_flux_y, first_flux_z, ok)
    if (.not. ok) return
    call distributed_fine_stage_3d( &
      distribution, old_state, old_state, first_flux_x, first_flux_y, &
      first_flux_z, dx, dy, dz, dt, .false., stage_state, ok)
    if (.not. ok) return
    allocate(stage_temperature(nx, ny, nz))
    call recover_reactive_temperatures_3d( &
      species, stage_state, old_temperature, nx, ny, nz, &
      stage_temperature, local_ok)
    call collective_logical_and(distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return

    allocate(second_flux_x(nvar, 0:nx, ny, nz))
    allocate(second_flux_y(nvar, nx, 0:ny, nz))
    allocate(second_flux_z(nvar, nx, ny, 0:nz))
    call distributed_fine_face_fluxes_3d( &
      distribution, species, patch, coarse_start, &
      coarse_start_temperature, coarse_end, coarse_end_temperature, &
      stage_state, stage_temperature, alpha_end, riemann_solver, &
      reconstruction, limiter, second_flux_x, second_flux_y, &
      second_flux_z, ok)
    if (.not. ok) return
    call distributed_fine_stage_3d( &
      distribution, old_state, stage_state, second_flux_x, second_flux_y, &
      second_flux_z, dx, dy, dz, dt, .true., candidate_state, ok)
    if (.not. ok) return
    allocate(candidate_temperature(nx, ny, nz))
    call recover_reactive_temperatures_3d( &
      species, candidate_state, stage_temperature, nx, ny, nz, &
      candidate_temperature, local_ok)
    call collective_logical_and(distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return

    face_flux_x = 0.5_dp * (first_flux_x + second_flux_x)
    face_flux_y = 0.5_dp * (first_flux_y + second_flux_y)
    face_flux_z = 0.5_dp * (first_flux_z + second_flux_z)
    state = candidate_state
    temperature = candidate_temperature
    ok = .true.
  end subroutine distributed_fine_patch_ssprk2_3d

  subroutine distributed_fine_face_fluxes_3d( &
      distribution, species, patch, coarse_start, &
      coarse_start_temperature, coarse_end, coarse_end_temperature, &
      state, temperature, alpha, riemann_solver, reconstruction, limiter, &
      flux_x, flux_y, flux_z, ok)
    type(mpi_amr_slab_distribution_3d), intent(in) :: distribution
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: coarse_start(:, :, :, :)
    real(dp), intent(in) :: coarse_start_temperature(:, :, :)
    real(dp), intent(in) :: coarse_end(:, :, :, :)
    real(dp), intent(in) :: coarse_end_temperature(:, :, :)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    real(dp), intent(in) :: alpha
    character(len=*), intent(in) :: riemann_solver, reconstruction, limiter
    real(dp), intent(out) :: flux_x(:, 0:, :, :)
    real(dp), intent(out) :: flux_y(:, :, 0:, :)
    real(dp), intent(out) :: flux_z(:, :, :, 0:)
    logical, intent(out) :: ok

    real(dp), allocatable :: local_x(:, :, :, :)
    real(dp), allocatable :: local_y(:, :, :, :)
    real(dp), allocatable :: local_z(:, :, :, :)
    real(dp), allocatable :: ghost_state(:), primitive(:)
    real(dp) :: ghost_temperature
    logical :: face_ok, local_ok, global_ok
    integer :: i, j, k, coarse_i, coarse_j, coarse_k
    integer :: nx, ny, nz, ratio

    nx = patch%fine_nx()
    ny = patch%fine_ny()
    nz = patch%fine_nz()
    ratio = patch%refinement_ratio
    allocate(local_x(size(state, 1), 0:nx, ny, nz))
    allocate(local_y(size(state, 1), nx, 0:ny, nz))
    allocate(local_z(size(state, 1), nx, ny, 0:nz))
    allocate(ghost_state(size(state, 1)))
    allocate(primitive(reactive_nprim(size(species))))
    local_x = 0.0_dp
    local_y = 0.0_dp
    local_z = 0.0_dp
    flux_x = 0.0_dp
    flux_y = 0.0_dp
    flux_z = 0.0_dp
    local_ok = ieee_is_finite(alpha)
    if (local_ok) local_ok = alpha >= 0.0_dp .and. alpha <= 1.0_dp
    if (local_ok) local_ok = &
      trim(reconstruction) == "pcm" .or. &
      trim(reconstruction) == "characteristic_plm"
    if (local_ok) local_ok = &
      trim(limiter) == "minmod" .or. trim(limiter) == "mc"
    if (local_ok .and. trim(reconstruction) == "characteristic_plm") then
      call compute_fine_patch_plm_slab_face_fluxes_3d( &
        species, patch, coarse_start, coarse_start_temperature, coarse_end, &
        coarse_end_temperature, state, temperature, alpha, &
        distribution%fine_first, distribution%fine_last, limiter, &
        riemann_solver, local_x, local_y, local_z, face_ok)
      local_ok = face_ok
    end if
    if (local_ok .and. trim(reconstruction) == "pcm") then
      x_faces: do k = 1, nz
        coarse_k = patch%coarse_k_lower + (k - 1) / ratio
        do j = 1, ny
          coarse_j = patch%coarse_j_lower + (j - 1) / ratio
          if (distribution%fine_first == 1) then
            call interpolate_distributed_coarse_cell( &
              species, coarse_start, coarse_start_temperature, &
              coarse_end, coarse_end_temperature, &
              patch%coarse_i_lower - 1, coarse_j, coarse_k, alpha, &
              primitive, ghost_state, ghost_temperature, face_ok)
            if (.not. face_ok) then
              local_ok = .false.
              exit x_faces
            end if
            call reactive_riemann_flux_x( &
              species, ghost_state, state(:, 1, j, k), &
              ghost_temperature, temperature(1, j, k), riemann_solver, &
              local_x(:, 0, j, k), face_ok)
            if (.not. face_ok) then
              local_ok = .false.
              exit x_faces
            end if
          end if
          do i = distribution%fine_first, distribution%fine_last
            if (i < nx) then
              call reactive_riemann_flux_x( &
                species, state(:, i, j, k), state(:, i + 1, j, k), &
                temperature(i, j, k), temperature(i + 1, j, k), &
                riemann_solver, local_x(:, i, j, k), face_ok)
            else
              call interpolate_distributed_coarse_cell( &
                species, coarse_start, coarse_start_temperature, &
                coarse_end, coarse_end_temperature, &
                patch%coarse_i_upper + 1, coarse_j, coarse_k, alpha, &
                primitive, ghost_state, ghost_temperature, face_ok)
              if (face_ok) then
                call reactive_riemann_flux_x( &
                  species, state(:, nx, j, k), ghost_state, &
                  temperature(nx, j, k), ghost_temperature, &
                  riemann_solver, local_x(:, nx, j, k), face_ok)
              end if
            end if
            if (.not. face_ok) then
              local_ok = .false.
              exit x_faces
            end if
          end do
        end do
      end do x_faces
    end if

    if (local_ok .and. trim(reconstruction) == "pcm") then
      y_faces: do k = 1, nz
        coarse_k = patch%coarse_k_lower + (k - 1) / ratio
        do i = distribution%fine_first, distribution%fine_last
          coarse_i = patch%coarse_i_lower + (i - 1) / ratio
          call interpolate_distributed_coarse_cell( &
            species, coarse_start, coarse_start_temperature, &
            coarse_end, coarse_end_temperature, coarse_i, &
            patch%coarse_j_lower - 1, coarse_k, alpha, primitive, &
            ghost_state, ghost_temperature, face_ok)
          if (.not. face_ok) then
            local_ok = .false.
            exit y_faces
          end if
          call reactive_riemann_flux_y( &
            species, ghost_state, state(:, i, 1, k), ghost_temperature, &
            temperature(i, 1, k), riemann_solver, &
            local_y(:, i, 0, k), face_ok)
          if (.not. face_ok) then
            local_ok = .false.
            exit y_faces
          end if
          do j = 1, ny - 1
            call reactive_riemann_flux_y( &
              species, state(:, i, j, k), state(:, i, j + 1, k), &
              temperature(i, j, k), temperature(i, j + 1, k), &
              riemann_solver, local_y(:, i, j, k), face_ok)
            if (.not. face_ok) then
              local_ok = .false.
              exit y_faces
            end if
          end do
          call interpolate_distributed_coarse_cell( &
            species, coarse_start, coarse_start_temperature, &
            coarse_end, coarse_end_temperature, coarse_i, &
            patch%coarse_j_upper + 1, coarse_k, alpha, primitive, &
            ghost_state, ghost_temperature, face_ok)
          if (face_ok) then
            call reactive_riemann_flux_y( &
              species, state(:, i, ny, k), ghost_state, &
              temperature(i, ny, k), ghost_temperature, riemann_solver, &
              local_y(:, i, ny, k), face_ok)
          end if
          if (.not. face_ok) then
            local_ok = .false.
            exit y_faces
          end if
        end do
      end do y_faces
    end if

    if (local_ok .and. trim(reconstruction) == "pcm") then
      z_faces: do j = 1, ny
        coarse_j = patch%coarse_j_lower + (j - 1) / ratio
        do i = distribution%fine_first, distribution%fine_last
          coarse_i = patch%coarse_i_lower + (i - 1) / ratio
          call interpolate_distributed_coarse_cell( &
            species, coarse_start, coarse_start_temperature, &
            coarse_end, coarse_end_temperature, coarse_i, coarse_j, &
            patch%coarse_k_lower - 1, alpha, primitive, ghost_state, &
            ghost_temperature, face_ok)
          if (.not. face_ok) then
            local_ok = .false.
            exit z_faces
          end if
          call reactive_riemann_flux_z( &
            species, ghost_state, state(:, i, j, 1), ghost_temperature, &
            temperature(i, j, 1), riemann_solver, &
            local_z(:, i, j, 0), face_ok)
          if (.not. face_ok) then
            local_ok = .false.
            exit z_faces
          end if
          do k = 1, nz - 1
            call reactive_riemann_flux_z( &
              species, state(:, i, j, k), state(:, i, j, k + 1), &
              temperature(i, j, k), temperature(i, j, k + 1), &
              riemann_solver, local_z(:, i, j, k), face_ok)
            if (.not. face_ok) then
              local_ok = .false.
              exit z_faces
            end if
          end do
          call interpolate_distributed_coarse_cell( &
            species, coarse_start, coarse_start_temperature, &
            coarse_end, coarse_end_temperature, coarse_i, coarse_j, &
            patch%coarse_k_upper + 1, alpha, primitive, ghost_state, &
            ghost_temperature, face_ok)
          if (face_ok) then
            call reactive_riemann_flux_z( &
              species, state(:, i, j, nz), ghost_state, &
              temperature(i, j, nz), ghost_temperature, riemann_solver, &
              local_z(:, i, j, nz), face_ok)
          end if
          if (.not. face_ok) then
            local_ok = .false.
            exit z_faces
          end if
        end do
      end do z_faces
    end if

    call collective_logical_and(distribution%comm, local_ok, global_ok)
    if (.not. global_ok) then
      ok = .false.
      return
    end if
    call collective_sum_rank4(distribution%comm, local_x, flux_x, ok)
    if (.not. ok) return
    call collective_sum_rank4(distribution%comm, local_y, flux_y, ok)
    if (.not. ok) return
    call collective_sum_rank4(distribution%comm, local_z, flux_z, ok)
  end subroutine distributed_fine_face_fluxes_3d

  subroutine interpolate_distributed_coarse_cell( &
      species, coarse_start, coarse_start_temperature, coarse_end, &
      coarse_end_temperature, i, j, k, alpha, primitive, state, &
      temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: coarse_start(:, :, :, :)
    real(dp), intent(in) :: coarse_start_temperature(:, :, :)
    real(dp), intent(in) :: coarse_end(:, :, :, :)
    real(dp), intent(in) :: coarse_end_temperature(:, :, :)
    integer, intent(in) :: i, j, k
    real(dp), intent(in) :: alpha
    real(dp), intent(out) :: primitive(:), state(:), temperature
    logical, intent(out) :: ok

    real(dp) :: temperature_guess, sound_speed

    state = (1.0_dp - alpha) * coarse_start(:, i, j, k) + &
      alpha * coarse_end(:, i, j, k)
    temperature_guess = &
      (1.0_dp - alpha) * coarse_start_temperature(i, j, k) + &
      alpha * coarse_end_temperature(i, j, k)
    call reactive_conserved_to_primitive( &
      species, state, temperature_guess, primitive, temperature, &
      sound_speed, ok)
  end subroutine interpolate_distributed_coarse_cell

  subroutine distributed_fine_stage_3d( &
      distribution, old_state, stage_base, flux_x, flux_y, flux_z, &
      dx, dy, dz, dt, second_stage, candidate, ok)
    type(mpi_amr_slab_distribution_3d), intent(in) :: distribution
    real(dp), intent(in) :: old_state(:, :, :, :)
    real(dp), intent(in) :: stage_base(:, :, :, :)
    real(dp), intent(in) :: flux_x(:, 0:, :, :)
    real(dp), intent(in) :: flux_y(:, :, 0:, :)
    real(dp), intent(in) :: flux_z(:, :, :, 0:)
    real(dp), intent(in) :: dx, dy, dz, dt
    logical, intent(in) :: second_stage
    real(dp), allocatable, intent(out) :: candidate(:, :, :, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: local_candidate(:, :, :, :), rhs(:)
    logical :: local_ok, global_ok
    integer :: i, j, k

    allocate(local_candidate, mold=old_state)
    allocate(candidate, mold=old_state)
    allocate(rhs(size(old_state, 1)))
    local_candidate = 0.0_dp
    candidate = 0.0_dp
    local_ok = distribution%fine_first >= 1 .and. &
      distribution%fine_last <= size(old_state, 2)
    if (local_ok) then
      do k = 1, size(old_state, 4)
        do j = 1, size(old_state, 3)
          do i = distribution%fine_first, distribution%fine_last
            rhs = -(flux_x(:, i, j, k) - flux_x(:, i - 1, j, k)) / dx &
              -(flux_y(:, i, j, k) - flux_y(:, i, j - 1, k)) / dy &
              -(flux_z(:, i, j, k) - flux_z(:, i, j, k - 1)) / dz
            if (second_stage) then
              local_candidate(:, i, j, k) = 0.5_dp * &
                old_state(:, i, j, k) + 0.5_dp * &
                  (stage_base(:, i, j, k) + dt * rhs)
            else
              local_candidate(:, i, j, k) = &
                stage_base(:, i, j, k) + dt * rhs
            end if
          end do
        end do
      end do
      local_ok = all(ieee_is_finite(local_candidate(:, &
        distribution%fine_first:distribution%fine_last, :, :)))
    end if
    call collective_logical_and(distribution%comm, local_ok, global_ok)
    if (.not. global_ok) then
      ok = .false.
      return
    end if
    call collective_sum_rank4( &
      distribution%comm, local_candidate, candidate, ok)
  end subroutine distributed_fine_stage_3d

  subroutine broadcast_mpi_amr_reactive_hierarchy_3d( &
      distribution, patch, root, coarse_state, coarse_temperature, &
      fine_state, fine_temperature, time, steps, initial_integrals, &
      maximum_reflux, ok)
    type(mpi_amr_slab_distribution_3d), intent(in) :: distribution
    type(amr_patch_3d), intent(in) :: patch
    integer, intent(in) :: root
    real(dp), intent(inout) :: coarse_state(:, :, :, :)
    real(dp), intent(inout) :: coarse_temperature(:, :, :)
    real(dp), intent(inout) :: fine_state(:, :, :, :)
    real(dp), intent(inout) :: fine_temperature(:, :, :)
    real(dp), intent(inout) :: time
    integer, intent(inout) :: steps
    real(dp), intent(inout) :: initial_integrals(:), maximum_reflux
    logical, intent(out) :: ok

    real(dp), allocatable :: candidate_coarse(:, :, :, :)
    real(dp), allocatable :: candidate_coarse_temperature(:, :, :)
    real(dp), allocatable :: candidate_fine(:, :, :, :)
    real(dp), allocatable :: candidate_fine_temperature(:, :, :)
    real(dp), allocatable :: candidate_integrals(:)
    real(dp) :: candidate_metadata(2)
    integer :: candidate_steps, ierr
    integer :: local_layout(25), root_layout(25)
    integer :: minimum_root, maximum_root
    logical :: local_ok, global_ok

    ok = .false.
    call MPI_Allreduce( &
      root, minimum_root, 1, MPI_INTEGER, MPI_MIN, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Allreduce( &
      root, maximum_root, 1, MPI_INTEGER, MPI_MAX, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    local_ok = distribution%is_valid(patch) .and. &
      root == minimum_root .and. root == maximum_root .and. &
      root >= 0 .and. root < distribution%nranks
    call collective_logical_and(distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    local_layout = [ &
      patch%coarse_nx, patch%coarse_ny, patch%coarse_nz, &
      patch%coarse_i_lower, patch%coarse_i_upper, &
      patch%coarse_j_lower, patch%coarse_j_upper, &
      patch%coarse_k_lower, patch%coarse_k_upper, patch%refinement_ratio, &
      shape(coarse_state), shape(coarse_temperature), shape(fine_state), &
      shape(fine_temperature), size(initial_integrals)]
    root_layout = 0
    if (distribution%rank == root) root_layout = local_layout
    call MPI_Bcast( &
      root_layout, size(root_layout), MPI_INTEGER, root, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    local_ok = all(local_layout == root_layout) .and. &
      size(initial_integrals) >= 1 .and. &
      all(shape(coarse_state) == [ &
        size(initial_integrals), patch%coarse_nx, patch%coarse_ny, &
        patch%coarse_nz]) .and. &
      all(shape(coarse_temperature) == [ &
        patch%coarse_nx, patch%coarse_ny, patch%coarse_nz]) .and. &
      all(shape(fine_state) == [ &
        size(initial_integrals), patch%fine_nx(), patch%fine_ny(), &
        patch%fine_nz()]) .and. &
      all(shape(fine_temperature) == [ &
        patch%fine_nx(), patch%fine_ny(), patch%fine_nz()])
    call collective_logical_and(distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    allocate(candidate_coarse, mold=coarse_state)
    allocate(candidate_coarse_temperature, mold=coarse_temperature)
    allocate(candidate_fine, mold=fine_state)
    allocate(candidate_fine_temperature, mold=fine_temperature)
    allocate(candidate_integrals, mold=initial_integrals)
    candidate_coarse = 0.0_dp
    candidate_coarse_temperature = 0.0_dp
    candidate_fine = 0.0_dp
    candidate_fine_temperature = 0.0_dp
    candidate_integrals = 0.0_dp
    candidate_metadata = 0.0_dp
    candidate_steps = 0
    if (distribution%rank == root) then
      candidate_coarse = coarse_state
      candidate_coarse_temperature = coarse_temperature
      candidate_fine = fine_state
      candidate_fine_temperature = fine_temperature
      candidate_integrals = initial_integrals
      candidate_metadata = [time, maximum_reflux]
      candidate_steps = steps
    end if
    call MPI_Bcast( &
      candidate_coarse, size(candidate_coarse), MPI_DOUBLE_PRECISION, &
      root, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Bcast( &
      candidate_coarse_temperature, size(candidate_coarse_temperature), &
      MPI_DOUBLE_PRECISION, root, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Bcast( &
      candidate_fine, size(candidate_fine), MPI_DOUBLE_PRECISION, &
      root, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Bcast( &
      candidate_fine_temperature, size(candidate_fine_temperature), &
      MPI_DOUBLE_PRECISION, root, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Bcast( &
      candidate_integrals, size(candidate_integrals), &
      MPI_DOUBLE_PRECISION, root, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Bcast( &
      candidate_metadata, 2, MPI_DOUBLE_PRECISION, root, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Bcast( &
      candidate_steps, 1, MPI_INTEGER, root, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) return
    local_ok = all(ieee_is_finite(candidate_coarse)) .and. &
      all(ieee_is_finite(candidate_coarse_temperature)) .and. &
      all(ieee_is_finite(candidate_fine)) .and. &
      all(ieee_is_finite(candidate_fine_temperature)) .and. &
      all(ieee_is_finite(candidate_integrals)) .and. &
      all(ieee_is_finite(candidate_metadata))
    if (local_ok) local_ok = &
      minval(candidate_coarse_temperature) > 0.0_dp .and. &
      minval(candidate_fine_temperature) > 0.0_dp .and. &
      candidate_metadata(1) >= 0.0_dp .and. &
      candidate_metadata(2) >= 0.0_dp .and. candidate_steps >= 0
    call collective_logical_and(distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    coarse_state = candidate_coarse
    coarse_temperature = candidate_coarse_temperature
    fine_state = candidate_fine
    fine_temperature = candidate_fine_temperature
    initial_integrals = candidate_integrals
    time = candidate_metadata(1)
    maximum_reflux = candidate_metadata(2)
    steps = candidate_steps
    ok = .true.
  end subroutine broadcast_mpi_amr_reactive_hierarchy_3d

  subroutine collective_metadata_matches( &
      distribution, species, patch, local_reals, local_solver, &
      local_valid, ok)
    type(mpi_amr_slab_distribution_3d), intent(in) :: distribution
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: local_reals(:)
    character(len=*), intent(in) :: local_solver
    logical, intent(in) :: local_valid
    logical, intent(out) :: ok

    real(dp), allocatable :: root_reals(:)
    real(dp) :: local_species_reals( &
      nasa7_real_count, maximum_species_count)
    real(dp) :: root_species_reals( &
      nasa7_real_count, maximum_species_count)
    integer :: local_patch(11), root_patch(11), ierr, species_index
    character(len=24) :: local_species_names(maximum_species_count)
    character(len=24) :: root_species_names(maximum_species_count)
    character(len=256) :: root_solver
    logical :: matches, global_matches

    allocate(root_reals(size(local_reals)))
    root_reals = local_reals
    local_patch = [ &
      patch%coarse_nx, patch%coarse_ny, patch%coarse_nz, &
      patch%coarse_i_lower, patch%coarse_i_upper, &
      patch%coarse_j_lower, patch%coarse_j_upper, &
      patch%coarse_k_lower, patch%coarse_k_upper, patch%refinement_ratio, &
      size(species)]
    root_patch = local_patch
    local_species_reals = 0.0_dp
    local_species_names = ""
    if (size(species) <= maximum_species_count) then
      do species_index = 1, size(species)
        local_species_names(species_index) = species(species_index)%name
        local_species_reals(:, species_index) = [ &
          species(species_index)%molecular_weight, &
          species(species_index)%temperature_min, &
          species(species_index)%temperature_mid, &
          species(species_index)%temperature_max, &
          species(species_index)%low_coefficients, &
          species(species_index)%high_coefficients]
      end do
    end if
    root_species_reals = local_species_reals
    root_species_names = local_species_names
    root_solver = ""
    if (distribution%rank == 0) root_solver = trim(local_solver)
    call MPI_Bcast( &
      root_reals, size(root_reals), MPI_DOUBLE_PRECISION, 0, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Bcast( &
      root_patch, size(root_patch), MPI_INTEGER, 0, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Bcast( &
      root_species_reals, size(root_species_reals), &
      MPI_DOUBLE_PRECISION, 0, distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Bcast( &
      root_species_names, len(root_species_names) * &
        size(root_species_names), MPI_CHARACTER, 0, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    call MPI_Bcast( &
      root_solver, len(root_solver), MPI_CHARACTER, 0, &
      distribution%comm, ierr)
    if (ierr /= MPI_SUCCESS) then
      ok = .false.
      return
    end if
    matches = local_valid .and. &
      all(identical_real_bits(local_reals, root_reals)) .and. &
      all(local_patch == root_patch) .and. &
      all(identical_real_bits( &
        local_species_reals, root_species_reals)) .and. &
      all(local_species_names == root_species_names) .and. &
      trim(local_solver) == trim(root_solver)
    call collective_logical_and( &
      distribution%comm, matches, global_matches)
    ok = global_matches
  end subroutine collective_metadata_matches

  subroutine collective_logical_and(comm, local_value, global_value)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(in) :: local_value
    logical, intent(out) :: global_value

    integer :: ierr

    call MPI_Allreduce( &
      local_value, global_value, 1, MPI_LOGICAL, MPI_LAND, comm, ierr)
    if (ierr /= MPI_SUCCESS) global_value = .false.
  end subroutine collective_logical_and

  pure elemental logical function identical_real_bits(left, right) &
      result(identical)
    real(dp), intent(in) :: left, right

    identical = transfer(left, 0_int64) == transfer(right, 0_int64)
  end function identical_real_bits

  subroutine collective_sum_rank4(comm, local_values, global_values, ok)
    type(MPI_Comm), intent(in) :: comm
    real(dp), intent(in), contiguous :: local_values(:, :, :, :)
    real(dp), intent(out), contiguous :: global_values(:, :, :, :)
    logical, intent(out) :: ok

    integer :: ierr

    call MPI_Allreduce( &
      local_values, global_values, size(local_values), &
      MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
    ok = ierr == MPI_SUCCESS
  end subroutine collective_sum_rank4

  pure logical function valid_mpi_amr_reactive_shapes( &
      species, patch, coarse_state, coarse_temperature, &
      fine_state, fine_temperature) result(valid)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: coarse_state(:, :, :, :)
    real(dp), intent(in) :: coarse_temperature(:, :, :)
    real(dp), intent(in) :: fine_state(:, :, :, :)
    real(dp), intent(in) :: fine_temperature(:, :, :)

    integer :: nvar

    nvar = reactive_nvar(size(species))
    valid = patch%is_strictly_interior() .and. &
      size(coarse_state, 1) == nvar .and. &
      all(shape(coarse_state) == [ &
        nvar, patch%coarse_nx, patch%coarse_ny, patch%coarse_nz]) .and. &
      all(shape(coarse_temperature) == [ &
        patch%coarse_nx, patch%coarse_ny, patch%coarse_nz]) .and. &
      all(shape(fine_state) == [ &
        nvar, patch%fine_nx(), patch%fine_ny(), patch%fine_nz()]) .and. &
      all(shape(fine_temperature) == [ &
        patch%fine_nx(), patch%fine_ny(), patch%fine_nz()]) .and. &
      all(ieee_is_finite(coarse_state)) .and. &
      all(ieee_is_finite(coarse_temperature)) .and. &
      all(ieee_is_finite(fine_state)) .and. &
      all(ieee_is_finite(fine_temperature))
  end function valid_mpi_amr_reactive_shapes

end module mpi_amr_reactive_3d_mod
