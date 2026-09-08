program test_mpi_amr_reactive_chemistry_3d
  use, intrinsic :: iso_fortran_env, only: int64
  use mpi_f08
  use precision_mod, only: dp
  use state_indices_mod, only: ncons
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: &
    load_h2o2_elementary_mechanism
  use mesh_3d_mod, only: uniform_cell_centers_3d
  use reactive_1d_mod, only: reactive_nvar
  use simulation_config_reactive_3d_mod, only: reactive_3d_config
  use reactive_entropy_wave_3d_problem_mod, only: &
    initialize_reactive_problem_3d
  use reactive_3d_mod, only: recover_reactive_temperatures_3d
  use amr_hierarchy_3d_mod, only: &
    amr_patch_3d, initialize_amr_patch_3d, average_down_3d
  use amr_reactive_3d_mod, only: &
    composite_element_integrals_amr_3d, &
    advance_amr_reactive_chemistry_3d, &
    advance_amr_reactive_strang_3d, advance_amr_reactive_hydro_3d
  use mpi_amr_sparse_reactive_3d_mod, only: &
    mpi_amr_sparse_distribution_3d, mpi_amr_sparse_hierarchy_3d, &
    initialize_mpi_amr_sparse_distribution_3d, &
    scatter_mpi_amr_sparse_hierarchy_3d, &
    gather_mpi_amr_sparse_hierarchy_3d, &
    advance_mpi_amr_sparse_chemistry_3d, &
    advance_mpi_amr_sparse_strang_3d
  implicit none

  integer, parameter :: nx = 8, ny = 4, nz = 4, ratio = 2
  real(dp), parameter :: chemistry_interval = 1.0e-6_dp
  real(dp), parameter :: dx = 1.0e-3_dp
  real(dp), parameter :: dy = 1.0e-3_dp
  real(dp), parameter :: dz = 1.0e-3_dp
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(elementary_reaction), allocatable :: mismatched_reactions(:)
  type(amr_patch_3d) :: patch
  type(mpi_amr_sparse_distribution_3d) :: distribution
  type(mpi_amr_sparse_hierarchy_3d) :: hierarchy, saved_hierarchy
  type(mpi_amr_sparse_hierarchy_3d) :: rejected_hierarchy
  type(mpi_amr_sparse_hierarchy_3d) :: saved_rejected_hierarchy
  real(dp), allocatable :: coarse_state(:, :, :, :)
  real(dp), allocatable :: fine_state(:, :, :, :)
  real(dp), allocatable :: coarse_temperature(:, :, :)
  real(dp), allocatable :: fine_temperature(:, :, :)
  real(dp), allocatable :: serial_coarse(:, :, :, :)
  real(dp), allocatable :: serial_fine(:, :, :, :)
  real(dp), allocatable :: serial_coarse_temperature(:, :, :)
  real(dp), allocatable :: serial_fine_temperature(:, :, :)
  real(dp), allocatable :: gathered_coarse(:, :, :, :)
  real(dp), allocatable :: gathered_fine(:, :, :, :)
  real(dp), allocatable :: gathered_coarse_temperature(:, :, :)
  real(dp), allocatable :: gathered_fine_temperature(:, :, :)
  real(dp) :: initial_elements(3), final_elements(3), element_error
  real(dp) :: chemistry_change, serial_reflux, sparse_reflux
  logical :: ok, local_condition, rank_chemistry_enabled
  integer :: ierr

  call MPI_Init(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Init failed"
  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "elementary thermodynamics load")
  call load_h2o2_elementary_mechanism(reactions, ok)
  call require(ok, "elementary mechanism load")
  call initialize_amr_patch_3d( &
    nx, ny, nz, 3, 6, 2, 3, 2, 3, ratio, patch, ok)
  call require(ok .and. patch%is_strictly_interior(), &
    "strictly interior sparse chemistry patch")
  call initialize_mpi_amr_sparse_distribution_3d( &
    patch, MPI_COMM_WORLD, distribution, ok)
  call require(ok, "sparse chemistry distribution")

  if (distribution%rank == 0) then
    call initialize_uniform_hierarchy( &
      species, patch, coarse_state, coarse_temperature, &
      fine_state, fine_temperature, ok)
  else
    allocate(coarse_state(0, 0, 0, 0), coarse_temperature(0, 0, 0))
    allocate(fine_state(0, 0, 0, 0), fine_temperature(0, 0, 0))
    ok = .true.
  end if
  call require(ok, "root uniform hierarchy initialization")
  if (distribution%rank == 0) then
    call composite_element_integrals_amr_3d( &
      species, patch, coarse_state, fine_state, dx, dy, dz, &
      initial_elements, ok)
  else
    ok = .true.
  end if
  call require(ok, "root initial element integrals")
  call reset_sparse_hierarchy(ok)
  call require(ok, "initial sparse hierarchy scatter")

  if (distribution%rank == 0) then
    serial_coarse = coarse_state
    serial_coarse_temperature = coarse_temperature
    serial_fine = fine_state
    serial_fine_temperature = fine_temperature
    call advance_amr_reactive_chemistry_3d( &
      species, reactions, patch, serial_coarse, serial_coarse_temperature, &
      serial_fine, serial_fine_temperature, chemistry_interval, &
      2.0e-7_dp, 1.0e-12_dp, ok)
  else
    ok = .true.
  end if
  call require(ok, "serial chemistry reference")
  call advance_mpi_amr_sparse_chemistry_3d( &
    distribution, species, reactions, patch, hierarchy, &
    chemistry_interval, 2.0e-7_dp, 1.0e-12_dp, ok)
  call require(ok, "rank-local sparse chemistry source")
  call gather_sparse_hierarchy(ok)
  call require(ok, "sparse chemistry gather")
  local_condition = .true.
  if (distribution%rank == 0) then
    local_condition = rank4_bits_match(gathered_coarse, serial_coarse) .and. &
      rank3_bits_match( &
        gathered_coarse_temperature, serial_coarse_temperature) .and. &
      rank4_bits_match(gathered_fine, serial_fine) .and. &
      rank3_bits_match(gathered_fine_temperature, serial_fine_temperature)
  end if
  call require(local_condition, "serial/sparse chemistry bit parity")
  if (distribution%rank == 0) then
    chemistry_change = max(maxval(abs( &
      serial_coarse(ncons + 1:, :, :, :) - &
      coarse_state(ncons + 1:, :, :, :))), maxval(abs( &
      serial_fine(ncons + 1:, :, :, :) - &
      fine_state(ncons + 1:, :, :, :))))
    call composite_element_integrals_amr_3d( &
      species, patch, serial_coarse, serial_fine, dx, dy, dz, &
      final_elements, ok)
    if (ok) element_error = maxval(abs(final_elements - initial_elements) / &
      max(1.0e-30_dp, abs(initial_elements)))
    local_condition = ok .and. chemistry_change > 1.0e-12_dp .and. &
      element_error <= 2.0e-10_dp
  else
    local_condition = .true.
  end if
  call require(local_condition, "sparse chemistry activity and elements")

  saved_hierarchy = hierarchy
  call advance_mpi_amr_sparse_chemistry_3d( &
    distribution, species, reactions, patch, hierarchy, &
    chemistry_interval, 2.0e-7_dp, 1.0e-12_dp, ok, "invalid")
  call require(.not. ok .and. &
    sparse_hierarchies_bits_match(hierarchy, saved_hierarchy), &
    "collective invalid-integrator rollback")

  if (distribution%nranks > 1) then
    mismatched_reactions = reactions
    if (distribution%rank == distribution%nranks - 1) then
      mismatched_reactions(1)%forward_rate%pre_exponential = nearest( &
        mismatched_reactions(1)%forward_rate%pre_exponential, 1.0_dp)
    end if
    call advance_mpi_amr_sparse_chemistry_3d( &
      distribution, species, mismatched_reactions, patch, hierarchy, &
      chemistry_interval, 2.0e-7_dp, 1.0e-12_dp, ok)
    call require(.not. ok .and. &
      sparse_hierarchies_bits_match(hierarchy, saved_hierarchy), &
      "rank-dependent reaction metadata rollback")
  end if

  if (any(distribution%fine_counts == 0)) then
    rejected_hierarchy = hierarchy
    if (distribution%fine_count == 0) &
      rejected_hierarchy%coarse_state(1, 1, 1, 1) = -1.0_dp
    saved_rejected_hierarchy = rejected_hierarchy
    call advance_mpi_amr_sparse_chemistry_3d( &
      distribution, species, reactions, patch, rejected_hierarchy, &
      chemistry_interval, 2.0e-7_dp, 1.0e-12_dp, ok)
    call require(.not. ok .and. sparse_hierarchies_bits_match( &
        rejected_hierarchy, saved_rejected_hierarchy), &
      "zero-fine rank collective source-failure rollback")
  end if

  call reset_sparse_hierarchy(ok)
  call require(ok, "split-parity sparse reset")
  if (distribution%rank == 0) then
    serial_coarse = coarse_state
    serial_coarse_temperature = coarse_temperature
    serial_fine = fine_state
    serial_fine_temperature = fine_temperature
    call advance_amr_reactive_strang_3d( &
      species, reactions, patch, serial_coarse, serial_coarse_temperature, &
      serial_fine, serial_fine_temperature, dx, dy, dz, &
      chemistry_interval, "rusanov", .true., 2.0e-7_dp, 1.0e-12_dp, &
      serial_reflux, ok, "characteristic_plm", "mc")
  else
    ok = .true.
  end if
  call require(ok, "serial reaction-hydro-reaction reference")
  call advance_mpi_amr_sparse_strang_3d( &
    distribution, species, reactions, patch, hierarchy, dx, dy, dz, &
    chemistry_interval, "rusanov", "characteristic_plm", "mc", .true., &
    2.0e-7_dp, 1.0e-12_dp, sparse_reflux, ok)
  call require(ok, "rank-local reaction-hydro-reaction split")
  call gather_sparse_hierarchy(ok)
  call require(ok, "rank-local split gather")
  local_condition = .true.
  if (distribution%rank == 0) then
    local_condition = rank4_bits_match(gathered_coarse, serial_coarse) .and. &
      rank3_bits_match( &
        gathered_coarse_temperature, serial_coarse_temperature) .and. &
      rank4_bits_match(gathered_fine, serial_fine) .and. &
      rank3_bits_match(gathered_fine_temperature, serial_fine_temperature) &
      .and. same_real_bits(sparse_reflux, serial_reflux)
  end if
  call require(local_condition, "serial/sparse split bit parity")

  saved_hierarchy = hierarchy
  if (distribution%nranks > 1) then
    rank_chemistry_enabled = distribution%rank /= distribution%nranks - 1
    call advance_mpi_amr_sparse_strang_3d( &
      distribution, species, reactions, patch, hierarchy, dx, dy, dz, &
      chemistry_interval, "rusanov", "characteristic_plm", "mc", &
      rank_chemistry_enabled, 2.0e-7_dp, 1.0e-12_dp, sparse_reflux, ok)
    call require(.not. ok .and. same_real_bits(sparse_reflux, 0.0_dp) .and. &
      sparse_hierarchies_bits_match(hierarchy, saved_hierarchy), &
      "rank-dependent chemistry-enable rollback")
  end if

  call reset_sparse_hierarchy(ok)
  call require(ok, "disabled-path sparse reset")
  if (distribution%rank == 0) then
    serial_coarse = coarse_state
    serial_coarse_temperature = coarse_temperature
    serial_fine = fine_state
    serial_fine_temperature = fine_temperature
    call advance_amr_reactive_hydro_3d( &
      species, patch, serial_coarse, serial_coarse_temperature, serial_fine, &
      serial_fine_temperature, dx, dy, dz, 1.0e-7_dp, "rusanov", &
      serial_reflux, ok, "characteristic_plm", "mc")
  else
    ok = .true.
  end if
  call require(ok, "chemistry-disabled serial hydro reference")
  call advance_mpi_amr_sparse_strang_3d( &
    distribution, species, reactions, patch, hierarchy, dx, dy, dz, &
    1.0e-7_dp, "rusanov", "characteristic_plm", "mc", .false., &
    -1.0_dp, -1.0_dp, sparse_reflux, ok)
  call require(ok, "chemistry-disabled sparse wrapper")
  call gather_sparse_hierarchy(ok)
  call require(ok, "chemistry-disabled sparse gather")
  local_condition = .true.
  if (distribution%rank == 0) then
    local_condition = rank4_bits_match(gathered_coarse, serial_coarse) .and. &
      rank3_bits_match( &
        gathered_coarse_temperature, serial_coarse_temperature) .and. &
      rank4_bits_match(gathered_fine, serial_fine) .and. &
      rank3_bits_match(gathered_fine_temperature, serial_fine_temperature) &
      .and. same_real_bits(sparse_reflux, serial_reflux)
  end if
  call require(local_condition, "chemistry-disabled sparse hydro bit parity")

  if (distribution%rank == 0) then
    write(*, '(2(a,es24.16))') &
      "chemistry change=", chemistry_change, ", element error=", element_error
    write(*, '(a)') "test_mpi_amr_reactive_chemistry_3d: PASS"
  end if
  call MPI_Finalize(ierr)
  if (ierr /= MPI_SUCCESS) error stop "MPI_Finalize failed"

contains

  subroutine reset_sparse_hierarchy(ok)
    logical, intent(out) :: ok

    call scatter_mpi_amr_sparse_hierarchy_3d( &
      distribution, species, patch, 0, coarse_state, coarse_temperature, &
      fine_state, fine_temperature, hierarchy, ok)
  end subroutine reset_sparse_hierarchy

  subroutine gather_sparse_hierarchy(ok)
    logical, intent(out) :: ok

    call gather_mpi_amr_sparse_hierarchy_3d( &
      distribution, patch, 0, hierarchy, gathered_coarse, &
      gathered_coarse_temperature, gathered_fine, &
      gathered_fine_temperature, ok)
  end subroutine gather_sparse_hierarchy

  subroutine initialize_uniform_hierarchy( &
      species, patch, coarse_state, coarse_temperature, &
      fine_state, fine_temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), allocatable, intent(out) :: coarse_state(:, :, :, :)
    real(dp), allocatable, intent(out) :: coarse_temperature(:, :, :)
    real(dp), allocatable, intent(out) :: fine_state(:, :, :, :)
    real(dp), allocatable, intent(out) :: fine_temperature(:, :, :)
    logical, intent(out) :: ok

    type(reactive_3d_config) :: config, fine_config
    real(dp), allocatable :: x(:), y(:), z(:), xf(:), yf(:), zf(:)
    real(dp), allocatable :: mass_fractions(:), recovered(:, :, :)
    real(dp) :: local_dx, local_dy, local_dz, density, fine_density
    logical :: local_ok
    integer :: nvar

    ok = .false.
    config = reactive_3d_config()
    config%nx = patch%coarse_nx
    config%ny = patch%coarse_ny
    config%nz = patch%coarse_nz
    config%x_upper = real(config%nx, dp) * dx
    config%y_upper = real(config%ny, dp) * dy
    config%z_upper = real(config%nz, dp) * dz
    config%problem = "uniform_reactor"
    config%thermo_model = "elementary"
    config%initial_temperature = 1200.0_dp
    config%initial_velocity_x = 0.0_dp
    config%initial_velocity_y = 0.0_dp
    config%initial_velocity_z = 0.0_dp
    nvar = reactive_nvar(size(species))
    allocate(coarse_state(nvar, config%nx, config%ny, config%nz))
    allocate(coarse_temperature(config%nx, config%ny, config%nz))
    allocate(fine_state(nvar, patch%fine_nx(), patch%fine_ny(), &
      patch%fine_nz()))
    allocate(fine_temperature( &
      patch%fine_nx(), patch%fine_ny(), patch%fine_nz()))
    allocate(x(config%nx), y(config%ny), z(config%nz))
    allocate(xf(patch%fine_nx()), yf(patch%fine_ny()), zf(patch%fine_nz()))
    allocate(mass_fractions(size(species)))
    allocate(recovered(config%nx, config%ny, config%nz))
    call uniform_cell_centers_3d( &
      config%nx, config%ny, config%nz, config%x_lower, config%x_upper, &
      config%y_lower, config%y_upper, config%z_lower, config%z_upper, &
      x, y, z, local_dx, local_dy, local_dz)
    call initialize_reactive_problem_3d( &
      species, config, x, y, z, coarse_state, coarse_temperature, &
      density, mass_fractions, local_ok)
    if (.not. local_ok) return
    fine_config = config
    fine_config%nx = patch%fine_nx()
    fine_config%ny = patch%fine_ny()
    fine_config%nz = patch%fine_nz()
    fine_config%x_lower = config%x_lower + &
      real(patch%coarse_i_lower - 1, dp) * dx
    fine_config%x_upper = config%x_lower + &
      real(patch%coarse_i_upper, dp) * dx
    fine_config%y_lower = config%y_lower + &
      real(patch%coarse_j_lower - 1, dp) * dy
    fine_config%y_upper = config%y_lower + &
      real(patch%coarse_j_upper, dp) * dy
    fine_config%z_lower = config%z_lower + &
      real(patch%coarse_k_lower - 1, dp) * dz
    fine_config%z_upper = config%z_lower + &
      real(patch%coarse_k_upper, dp) * dz
    call uniform_cell_centers_3d( &
      fine_config%nx, fine_config%ny, fine_config%nz, &
      fine_config%x_lower, fine_config%x_upper, fine_config%y_lower, &
      fine_config%y_upper, fine_config%z_lower, fine_config%z_upper, &
      xf, yf, zf, local_dx, local_dy, local_dz)
    call initialize_reactive_problem_3d( &
      species, fine_config, xf, yf, zf, fine_state, fine_temperature, &
      fine_density, mass_fractions, local_ok)
    if (.not. local_ok) return
    call average_down_3d(coarse_state, fine_state, patch, local_ok)
    if (.not. local_ok) return
    call recover_reactive_temperatures_3d( &
      species, coarse_state, coarse_temperature, config%nx, config%ny, &
      config%nz, recovered, local_ok)
    if (.not. local_ok) return
    coarse_temperature = recovered
    ok = .true.
  end subroutine initialize_uniform_hierarchy

  pure logical function same_real_bits(left, right) result(matches)
    real(dp), intent(in) :: left, right

    matches = transfer(left, 0_int64) == transfer(right, 0_int64)
  end function same_real_bits

  pure logical function rank4_bits_match(left, right) result(matches)
    real(dp), intent(in) :: left(:, :, :, :), right(:, :, :, :)
    integer :: i, j, k, component

    matches = all(shape(left) == shape(right))
    if (.not. matches) return
    do k = 1, size(left, 4)
      do j = 1, size(left, 3)
        do i = 1, size(left, 2)
          do component = 1, size(left, 1)
            if (.not. same_real_bits( &
                left(component, i, j, k), right(component, i, j, k))) then
              matches = .false.
              return
            end if
          end do
        end do
      end do
    end do
  end function rank4_bits_match

  pure logical function rank3_bits_match(left, right) result(matches)
    real(dp), intent(in) :: left(:, :, :), right(:, :, :)
    integer :: i, j, k

    matches = all(shape(left) == shape(right))
    if (.not. matches) return
    do k = 1, size(left, 3)
      do j = 1, size(left, 2)
        do i = 1, size(left, 1)
          if (.not. same_real_bits(left(i, j, k), right(i, j, k))) then
            matches = .false.
            return
          end if
        end do
      end do
    end do
  end function rank3_bits_match

  pure logical function sparse_hierarchies_bits_match( &
      left, right) result(matches)
    type(mpi_amr_sparse_hierarchy_3d), intent(in) :: left, right

    matches = left%nvar == right%nvar .and. &
      rank4_bits_match(left%coarse_state, right%coarse_state) .and. &
      rank3_bits_match( &
        left%coarse_temperature, right%coarse_temperature) .and. &
      rank4_bits_match(left%fine_state, right%fine_state) .and. &
      rank3_bits_match(left%fine_temperature, right%fine_temperature)
  end function sparse_hierarchies_bits_match

  subroutine require(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label

    logical :: global_condition
    integer :: local_ierr

    call MPI_Allreduce(condition, global_condition, 1, MPI_LOGICAL, &
      MPI_LAND, MPI_COMM_WORLD, local_ierr)
    if (local_ierr /= MPI_SUCCESS .or. .not. global_condition) then
      if (distribution%rank <= 0) write(*, '(a,1x,a)') "FAILED:", trim(label)
      call MPI_Abort(MPI_COMM_WORLD, 1, local_ierr)
      error stop 1
    end if
  end subroutine require

end program test_mpi_amr_reactive_chemistry_3d
