module mpi_amr_sparse_reactive_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: int64
  use mpi_f08
  use precision_mod, only: dp
  use constants_mod, only: density_floor, pressure_floor
  use nasa7_thermo_mod, only: nasa7_species
  use gas_transport_mod, only: &
    gas_transport_species, compatible_transport_database
  use elementary_kinetics_mod, only: &
    elementary_reaction, valid_elementary_reaction
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_conserved_to_primitive, &
    reactive_primitive_to_conserved, reactive_mass_fraction_component, &
    reactive_riemann_flux_x
  use slope_limiter_mod, only: limited_slope
  use reactive_2d_mod, only: reactive_riemann_flux_y
  use reactive_directional_flux_3d_mod, only: reactive_riemann_flux_z
  use reactive_3d_mod, only: &
    recover_reactive_temperatures_3d, &
    compute_reactive_plm_slab_face_fluxes_3d, &
    advance_reactive_chemistry_3d
  use reactive_transport_3d_mod, only: &
    reactive_transport_ghosted_fluxes_3d, &
    reactive_transport_interface_theta_3d, reactive_transport_timestep_3d
  use amr_hierarchy_3d_mod, only: amr_patch_3d
  implicit none
  private

  integer, parameter, public :: mpi_amr_sparse_ghost_width_3d = 2
  integer, parameter :: sparse_halo_left_tag_3d = 3410
  integer, parameter :: sparse_halo_right_tag_3d = 3420
  integer, parameter :: sparse_theta_left_tag_3d = 3430
  integer, parameter :: sparse_theta_right_tag_3d = 3440
  integer, parameter :: sparse_maximum_species_count_3d = 32
  integer, parameter :: sparse_nasa7_real_count_3d = 18

  type, public :: mpi_amr_sparse_distribution_3d
    type(MPI_Comm) :: comm = MPI_COMM_NULL
    integer :: rank = -1
    integer :: nranks = 0
    integer :: coarse_first = 0
    integer :: coarse_last = -1
    integer :: coarse_count = 0
    integer :: fine_first = 0
    integer :: fine_last = -1
    integer :: fine_count = 0
    integer :: previous_fine_rank = MPI_PROC_NULL
    integer :: next_fine_rank = MPI_PROC_NULL
    integer, allocatable :: coarse_counts(:)
    integer, allocatable :: coarse_displacements(:)
    integer, allocatable :: fine_counts(:)
    integer, allocatable :: fine_displacements(:)
  contains
    procedure :: is_valid => sparse_distribution_is_valid_3d
    procedure :: owns_coarse => sparse_distribution_owns_coarse_3d
    procedure :: owns_fine => sparse_distribution_owns_fine_3d
    procedure :: coarse_local_index => sparse_coarse_local_index_3d
    procedure :: fine_local_index => sparse_fine_local_index_3d
  end type mpi_amr_sparse_distribution_3d

  type, public :: mpi_amr_sparse_hierarchy_3d
    integer :: nvar = 0
    real(dp), allocatable :: coarse_state(:, :, :, :)
    real(dp), allocatable :: coarse_temperature(:, :, :)
    real(dp), allocatable :: fine_state(:, :, :, :)
    real(dp), allocatable :: fine_temperature(:, :, :)
  contains
    procedure :: is_valid => sparse_hierarchy_is_valid_3d
    procedure :: local_cell_count => sparse_hierarchy_local_cell_count_3d
    procedure :: local_value_count => sparse_hierarchy_local_value_count_3d
  end type mpi_amr_sparse_hierarchy_3d

  public :: initialize_mpi_amr_sparse_distribution_3d
  public :: scatter_mpi_amr_sparse_hierarchy_3d
  public :: gather_mpi_amr_sparse_hierarchy_3d
  public :: compute_mpi_amr_sparse_cfl_timestep_3d
  public :: compute_mpi_amr_sparse_transport_timestep_3d
  public :: build_mpi_amr_sparse_periodic_halos_3d
  public :: advance_mpi_amr_sparse_periodic_hydro_3d
  public :: advance_mpi_amr_sparse_hydro_3d
  public :: advance_mpi_amr_sparse_chemistry_3d
  public :: advance_mpi_amr_sparse_strang_3d
  public :: advance_mpi_amr_sparse_transport_3d
  public :: advance_mpi_amr_sparse_full_3d

contains

  subroutine advance_mpi_amr_sparse_chemistry_3d( &
      distribution, species, reactions, patch, hierarchy, interval, &
      rtol, atol, ok, chemistry_integrator)
    type(mpi_amr_sparse_distribution_3d), intent(in) :: distribution
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(amr_patch_3d), intent(in) :: patch
    type(mpi_amr_sparse_hierarchy_3d), intent(inout) :: hierarchy
    real(dp), intent(in) :: interval, rtol, atol
    logical, intent(out) :: ok
    character(len=*), intent(in), optional :: chemistry_integrator

    type(mpi_amr_sparse_hierarchy_3d) :: candidate
    real(dp), allocatable :: synchronized_temperature(:, :, :)
    logical :: coarse_ok, fine_ok, local_ok, global_ok
    character(len=32) :: integrator_name

    ok = .false.
    integrator_name = "default"
    if (present(chemistry_integrator)) &
      integrator_name = trim(chemistry_integrator)
    local_ok = all(ieee_is_finite([interval, rtol, atol]))
    if (local_ok) local_ok = interval >= 0.0_dp .and. &
      rtol > 0.0_dp .and. atol > 0.0_dp
    if (local_ok) local_ok = hierarchy%nvar == reactive_nvar(size(species))
    if (local_ok) local_ok = hierarchy%is_valid(distribution, patch)
    if (local_ok) local_ok = size(reactions) >= 1
    call sparse_collective_contract_matches_3d( &
      distribution, species, patch, [interval, rtol, atol], &
      "chemistry|" // trim(integrator_name), local_ok, global_ok)
    if (.not. global_ok) return
    call sparse_collective_reactions_match_3d( &
      distribution, reactions, size(species), .true., global_ok)
    if (.not. global_ok) return

    candidate = hierarchy
    call advance_reactive_chemistry_3d( &
      species, reactions, candidate%coarse_state, &
      candidate%coarse_temperature, distribution%coarse_count, &
      patch%coarse_ny, patch%coarse_nz, interval, rtol, atol, coarse_ok, &
      chemistry_integrator=chemistry_integrator)
    if (distribution%fine_count > 0) then
      call advance_reactive_chemistry_3d( &
        species, reactions, candidate%fine_state, &
        candidate%fine_temperature, distribution%fine_count, &
        patch%fine_ny(), patch%fine_nz(), interval, rtol, atol, fine_ok, &
        chemistry_integrator=chemistry_integrator)
    else
      fine_ok = .true.
    end if
    local_ok = coarse_ok .and. fine_ok
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return

    call average_down_sparse_hierarchy_3d( &
      distribution, patch, candidate, local_ok)
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    allocate(synchronized_temperature, mold=candidate%coarse_temperature)
    call recover_reactive_temperatures_3d( &
      species, candidate%coarse_state, candidate%coarse_temperature, &
      distribution%coarse_count, patch%coarse_ny, patch%coarse_nz, &
      synchronized_temperature, local_ok)
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    candidate%coarse_temperature = synchronized_temperature
    local_ok = candidate%is_valid(distribution, patch)
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return

    hierarchy = candidate
    ok = .true.
  end subroutine advance_mpi_amr_sparse_chemistry_3d

  subroutine advance_mpi_amr_sparse_strang_3d( &
      distribution, species, reactions, patch, hierarchy, dx, dy, dz, dt, &
      riemann_solver, reconstruction, limiter, chemistry_enabled, rtol, &
      atol, maximum_reflux_correction, ok, chemistry_integrator)
    type(mpi_amr_sparse_distribution_3d), intent(in) :: distribution
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(amr_patch_3d), intent(in) :: patch
    type(mpi_amr_sparse_hierarchy_3d), intent(inout) :: hierarchy
    real(dp), intent(in) :: dx, dy, dz, dt, rtol, atol
    character(len=*), intent(in) :: riemann_solver, reconstruction, limiter
    logical, intent(in) :: chemistry_enabled
    real(dp), intent(out) :: maximum_reflux_correction
    logical, intent(out) :: ok
    character(len=*), intent(in), optional :: chemistry_integrator

    type(mpi_amr_sparse_hierarchy_3d) :: candidate
    real(dp) :: contract_rtol, contract_atol
    logical :: local_ok, global_ok
    character(len=32) :: integrator_name

    maximum_reflux_correction = 0.0_dp
    ok = .false.
    integrator_name = "disabled"
    contract_rtol = 0.0_dp
    contract_atol = 0.0_dp
    if (chemistry_enabled) then
      integrator_name = "default"
      if (present(chemistry_integrator)) &
        integrator_name = trim(chemistry_integrator)
      contract_rtol = rtol
      contract_atol = atol
    end if
    local_ok = all(ieee_is_finite([dx, dy, dz, dt]))
    if (local_ok) local_ok = dx > 0.0_dp .and. dy > 0.0_dp .and. &
      dz > 0.0_dp .and. dt > 0.0_dp
    if (local_ok) local_ok = hierarchy%nvar == reactive_nvar(size(species))
    if (local_ok) local_ok = hierarchy%is_valid(distribution, patch)
    if (local_ok) local_ok = trim(reconstruction) == "pcm" .or. &
      trim(reconstruction) == "characteristic_plm"
    if (local_ok) local_ok = trim(limiter) == "minmod" .or. &
      trim(limiter) == "mc"
    if (chemistry_enabled) then
      if (local_ok) local_ok = all(ieee_is_finite([rtol, atol]))
      if (local_ok) local_ok = size(reactions) >= 1
      if (local_ok) local_ok = rtol > 0.0_dp .and. atol > 0.0_dp
    end if
    call sparse_collective_contract_matches_3d( &
      distribution, species, patch, &
      [dx, dy, dz, dt, contract_rtol, contract_atol, &
        merge(1.0_dp, 0.0_dp, chemistry_enabled)], &
      trim(riemann_solver) // "|" // trim(reconstruction) // "|" // &
        trim(limiter) // "|" // trim(integrator_name), &
      local_ok, global_ok)
    if (.not. global_ok) return

    candidate = hierarchy
    if (chemistry_enabled) then
      call advance_mpi_amr_sparse_chemistry_3d( &
        distribution, species, reactions, patch, candidate, &
        0.5_dp * dt, rtol, atol, local_ok, chemistry_integrator)
      if (.not. local_ok) return
    end if
    call advance_mpi_amr_sparse_hydro_3d( &
      distribution, species, patch, candidate, dx, dy, dz, dt, &
      riemann_solver, reconstruction, limiter, &
      maximum_reflux_correction, local_ok)
    if (.not. local_ok) then
      maximum_reflux_correction = 0.0_dp
      return
    end if
    if (chemistry_enabled) then
      call advance_mpi_amr_sparse_chemistry_3d( &
        distribution, species, reactions, patch, candidate, &
        0.5_dp * dt, rtol, atol, local_ok, chemistry_integrator)
      if (.not. local_ok) then
        maximum_reflux_correction = 0.0_dp
        return
      end if
    end if

    hierarchy = candidate
    ok = .true.
  end subroutine advance_mpi_amr_sparse_strang_3d

  subroutine advance_mpi_amr_sparse_transport_3d( &
      distribution, species, transport, patch, hierarchy, dx, dy, dz, &
      interval, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, limiter, &
      minimum_theta, maximum_reflux_correction, ok)
    type(mpi_amr_sparse_distribution_3d), intent(in) :: distribution
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(amr_patch_3d), intent(in) :: patch
    type(mpi_amr_sparse_hierarchy_3d), intent(inout) :: hierarchy
    real(dp), intent(in) :: dx, dy, dz, interval
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    character(len=*), intent(in) :: limiter
    real(dp), intent(out) :: minimum_theta, maximum_reflux_correction
    logical, intent(out) :: ok

    type(mpi_amr_sparse_hierarchy_3d) :: stage, euler, candidate
    real(dp), allocatable :: recovered_temperature(:, :, :)
    real(dp) :: theta_one, theta_two, reflux_one, reflux_two
    logical :: local_ok, global_ok

    minimum_theta = 1.0_dp
    maximum_reflux_correction = 0.0_dp
    ok = .false.
    local_ok = all(ieee_is_finite([dx, dy, dz, interval]))
    if (local_ok) local_ok = dx > 0.0_dp .and. dy > 0.0_dp .and. &
      dz > 0.0_dp .and. interval >= 0.0_dp
    if (local_ok) local_ok = hierarchy%nvar == reactive_nvar(size(species))
    if (local_ok) local_ok = hierarchy%is_valid(distribution, patch)
    if (local_ok) local_ok = compatible_transport_database(species, transport)
    if (local_ok) local_ok = .not. barodiffusion_enabled .or. &
      species_diffusion_enabled
    if (local_ok) local_ok = trim(limiter) == "minmod" .or. &
      trim(limiter) == "mc"
    call sparse_collective_contract_matches_3d( &
      distribution, species, patch, &
      [dx, dy, dz, interval, &
        merge(1.0_dp, 0.0_dp, viscosity_enabled), &
        merge(1.0_dp, 0.0_dp, thermal_conduction_enabled), &
        merge(1.0_dp, 0.0_dp, species_diffusion_enabled), &
        merge(1.0_dp, 0.0_dp, barodiffusion_enabled)], &
      "transport|" // trim(limiter), local_ok, global_ok)
    if (.not. global_ok) return
    call sparse_collective_transport_matches_3d( &
      distribution, transport, .true., global_ok)
    if (.not. global_ok) return
    if (interval <= tiny(1.0_dp) .or. .not. (viscosity_enabled .or. &
        thermal_conduction_enabled .or. species_diffusion_enabled)) then
      ok = .true.
      return
    end if

    call advance_mpi_amr_sparse_transport_euler_3d( &
      distribution, species, transport, patch, hierarchy, dx, dy, dz, &
      interval, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, limiter, stage, &
      theta_one, reflux_one, local_ok)
    if (.not. local_ok) return
    call advance_mpi_amr_sparse_transport_euler_3d( &
      distribution, species, transport, patch, stage, dx, dy, dz, &
      interval, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, limiter, euler, &
      theta_two, reflux_two, local_ok)
    if (.not. local_ok) return

    candidate = hierarchy
    candidate%coarse_state = &
      0.5_dp * (hierarchy%coarse_state + euler%coarse_state)
    candidate%fine_state = &
      0.5_dp * (hierarchy%fine_state + euler%fine_state)
    if (distribution%fine_count > 0) then
      allocate(recovered_temperature, mold=candidate%fine_temperature)
      call recover_reactive_temperatures_3d( &
        species, candidate%fine_state, &
        0.5_dp * (hierarchy%fine_temperature + euler%fine_temperature), &
        distribution%fine_count, patch%fine_ny(), patch%fine_nz(), &
        recovered_temperature, local_ok)
      if (local_ok) candidate%fine_temperature = recovered_temperature
      deallocate(recovered_temperature)
    else
      local_ok = .true.
    end if
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    call average_down_sparse_hierarchy_3d( &
      distribution, patch, candidate, local_ok)
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    allocate(recovered_temperature, mold=candidate%coarse_temperature)
    call recover_reactive_temperatures_3d( &
      species, candidate%coarse_state, &
      0.5_dp * (hierarchy%coarse_temperature + euler%coarse_temperature), &
      distribution%coarse_count, patch%coarse_ny, patch%coarse_nz, &
      recovered_temperature, local_ok)
    if (local_ok) candidate%coarse_temperature = recovered_temperature
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    local_ok = candidate%is_valid(distribution, patch)
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return

    hierarchy = candidate
    minimum_theta = min(theta_one, theta_two)
    maximum_reflux_correction = max(reflux_one, reflux_two)
    ok = .true.
  end subroutine advance_mpi_amr_sparse_transport_3d

  subroutine advance_mpi_amr_sparse_full_3d( &
      distribution, species, reactions, transport, patch, hierarchy, &
      dx, dy, dz, dt, riemann_solver, reconstruction, limiter, &
      chemistry_enabled, rtol, atol, transport_enabled, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, minimum_transport_theta, &
      maximum_reflux_correction, ok, chemistry_integrator)
    type(mpi_amr_sparse_distribution_3d), intent(in) :: distribution
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(amr_patch_3d), intent(in) :: patch
    type(mpi_amr_sparse_hierarchy_3d), intent(inout) :: hierarchy
    real(dp), intent(in) :: dx, dy, dz, dt, rtol, atol
    character(len=*), intent(in) :: riemann_solver, reconstruction, limiter
    logical, intent(in) :: chemistry_enabled, transport_enabled
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), intent(out) :: minimum_transport_theta
    real(dp), intent(out) :: maximum_reflux_correction
    logical, intent(out) :: ok
    character(len=*), intent(in), optional :: chemistry_integrator

    type(mpi_amr_sparse_hierarchy_3d) :: candidate
    real(dp) :: stage_theta, stage_reflux, hydro_reflux
    real(dp) :: candidate_theta, candidate_reflux
    real(dp) :: contract_rtol, contract_atol
    logical :: local_ok, global_ok
    character(len=32) :: integrator_name

    minimum_transport_theta = 1.0_dp
    maximum_reflux_correction = 0.0_dp
    ok = .false.
    candidate_theta = 1.0_dp
    candidate_reflux = 0.0_dp
    integrator_name = "disabled"
    contract_rtol = 0.0_dp
    contract_atol = 0.0_dp
    if (chemistry_enabled) then
      integrator_name = "default"
      if (present(chemistry_integrator)) &
        integrator_name = trim(chemistry_integrator)
      contract_rtol = rtol
      contract_atol = atol
    end if
    local_ok = all(ieee_is_finite([dx, dy, dz, dt]))
    if (local_ok) local_ok = dx > 0.0_dp .and. dy > 0.0_dp .and. &
      dz > 0.0_dp .and. dt > 0.0_dp
    if (local_ok) local_ok = hierarchy%nvar == reactive_nvar(size(species))
    if (local_ok) local_ok = hierarchy%is_valid(distribution, patch)
    if (local_ok) local_ok = trim(reconstruction) == "pcm" .or. &
      trim(reconstruction) == "characteristic_plm"
    if (local_ok) local_ok = trim(limiter) == "minmod" .or. &
      trim(limiter) == "mc"
    if (chemistry_enabled) then
      if (local_ok) local_ok = all(ieee_is_finite([rtol, atol]))
      if (local_ok) local_ok = size(reactions) >= 1
      if (local_ok) local_ok = rtol > 0.0_dp .and. atol > 0.0_dp
    end if
    if (transport_enabled) then
      if (local_ok) local_ok = compatible_transport_database(species, transport)
      if (local_ok) local_ok = .not. barodiffusion_enabled .or. &
        species_diffusion_enabled
    end if
    call sparse_collective_contract_matches_3d( &
      distribution, species, patch, &
      [dx, dy, dz, dt, contract_rtol, contract_atol, &
        merge(1.0_dp, 0.0_dp, chemistry_enabled), &
        merge(1.0_dp, 0.0_dp, transport_enabled), &
        merge(1.0_dp, 0.0_dp, viscosity_enabled), &
        merge(1.0_dp, 0.0_dp, thermal_conduction_enabled), &
        merge(1.0_dp, 0.0_dp, species_diffusion_enabled), &
        merge(1.0_dp, 0.0_dp, barodiffusion_enabled)], &
      trim(riemann_solver) // "|" // trim(reconstruction) // "|" // &
        trim(limiter) // "|" // trim(integrator_name), &
      local_ok, global_ok)
    if (.not. global_ok) return
    if (.not. transport_enabled) then
      call advance_mpi_amr_sparse_strang_3d( &
        distribution, species, reactions, patch, hierarchy, dx, dy, dz, &
        dt, riemann_solver, reconstruction, limiter, chemistry_enabled, &
        rtol, atol, candidate_reflux, ok, chemistry_integrator)
      if (ok) maximum_reflux_correction = candidate_reflux
      return
    end if
    call sparse_collective_transport_matches_3d( &
      distribution, transport, .true., global_ok)
    if (.not. global_ok) return

    candidate = hierarchy
    if (chemistry_enabled) then
      call advance_mpi_amr_sparse_chemistry_3d( &
        distribution, species, reactions, patch, candidate, &
        0.5_dp * dt, rtol, atol, local_ok, chemistry_integrator)
      if (.not. local_ok) return
    end if
    call advance_mpi_amr_sparse_transport_3d( &
      distribution, species, transport, patch, candidate, dx, dy, dz, &
      0.5_dp * dt, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, limiter, &
      stage_theta, stage_reflux, local_ok)
    if (.not. local_ok) return
    candidate_theta = min(candidate_theta, stage_theta)
    candidate_reflux = max(candidate_reflux, stage_reflux)
    call advance_mpi_amr_sparse_hydro_3d( &
      distribution, species, patch, candidate, dx, dy, dz, dt, &
      riemann_solver, reconstruction, limiter, hydro_reflux, local_ok)
    if (.not. local_ok) return
    candidate_reflux = max(candidate_reflux, hydro_reflux)
    call advance_mpi_amr_sparse_transport_3d( &
      distribution, species, transport, patch, candidate, dx, dy, dz, &
      0.5_dp * dt, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, limiter, &
      stage_theta, stage_reflux, local_ok)
    if (.not. local_ok) return
    candidate_theta = min(candidate_theta, stage_theta)
    candidate_reflux = max(candidate_reflux, stage_reflux)
    if (chemistry_enabled) then
      call advance_mpi_amr_sparse_chemistry_3d( &
        distribution, species, reactions, patch, candidate, &
        0.5_dp * dt, rtol, atol, local_ok, chemistry_integrator)
      if (.not. local_ok) return
    end if

    hierarchy = candidate
    minimum_transport_theta = candidate_theta
    maximum_reflux_correction = candidate_reflux
    ok = .true.
  end subroutine advance_mpi_amr_sparse_full_3d

  subroutine advance_mpi_amr_sparse_transport_euler_3d( &
      distribution, species, transport, patch, hierarchy, dx, dy, dz, dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, limiter, &
      new_hierarchy, minimum_theta, maximum_reflux_correction, ok)
    type(mpi_amr_sparse_distribution_3d), intent(in) :: distribution
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(amr_patch_3d), intent(in) :: patch
    type(mpi_amr_sparse_hierarchy_3d), intent(in) :: hierarchy
    real(dp), intent(in) :: dx, dy, dz, dt
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    character(len=*), intent(in) :: limiter
    type(mpi_amr_sparse_hierarchy_3d), intent(out) :: new_hierarchy
    real(dp), intent(out) :: minimum_theta, maximum_reflux_correction
    logical, intent(out) :: ok

    real(dp), allocatable :: coarse_start_halo(:, :, :, :)
    real(dp), allocatable :: coarse_start_halo_temperature(:, :, :)
    real(dp), allocatable :: coarse_end_halo(:, :, :, :)
    real(dp), allocatable :: coarse_end_halo_temperature(:, :, :)
    real(dp), allocatable :: recovered_temperature(:, :, :)
    real(dp), allocatable :: coarse_flux_x(:, :, :, :)
    real(dp), allocatable :: coarse_flux_y(:, :, :, :)
    real(dp), allocatable :: coarse_flux_z(:, :, :, :)
    real(dp), allocatable :: fine_flux_x(:, :, :, :)
    real(dp), allocatable :: fine_flux_y(:, :, :, :)
    real(dp), allocatable :: fine_flux_z(:, :, :, :)
    real(dp), allocatable :: fine_x_lower(:, :, :), fine_x_upper(:, :, :)
    real(dp), allocatable :: fine_y_lower(:, :, :), fine_y_upper(:, :, :)
    real(dp), allocatable :: fine_z_lower(:, :, :), fine_z_upper(:, :, :)
    real(dp), allocatable :: boundary_x_lower(:, :, :)
    real(dp), allocatable :: boundary_x_upper(:, :, :)
    real(dp), allocatable :: boundary_y_lower(:, :, :)
    real(dp), allocatable :: boundary_y_upper(:, :, :)
    real(dp), allocatable :: boundary_z_lower(:, :, :)
    real(dp), allocatable :: boundary_z_upper(:, :, :)
    real(dp) :: alpha, fine_dt, ratio_real
    real(dp) :: coarse_theta, fine_theta, local_reflux, global_reflux
    real(dp) :: interface_theta
    logical :: local_ok, global_ok
    integer :: i, j, k, previous_j, previous_k
    integer :: ratio, subcycles, substep, covered_nx, covered_ny, covered_nz
    integer :: ierr

    new_hierarchy = hierarchy
    minimum_theta = 1.0_dp
    maximum_reflux_correction = 0.0_dp
    ok = .false.

    call compute_sparse_periodic_transport_fluxes_3d( &
      distribution, species, transport, patch, hierarchy%coarse_state, &
      hierarchy%coarse_temperature, dx, dy, dz, dt, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, coarse_flux_x, coarse_flux_y, coarse_flux_z, &
      coarse_theta, local_ok)
    if (.not. local_ok) return
    do k = 1, patch%coarse_nz
      previous_k = 1 + modulo(k - 2, patch%coarse_nz)
      do j = 1, patch%coarse_ny
        previous_j = 1 + modulo(j - 2, patch%coarse_ny)
        do i = 1, distribution%coarse_count
          new_hierarchy%coarse_state(:, i, j, k) = &
            hierarchy%coarse_state(:, i, j, k) - &
            dt / dx * (coarse_flux_x(:, i, j, k) - &
              coarse_flux_x(:, i - 1, j, k)) - &
            dt / dy * (coarse_flux_y(:, i, j, k) - &
              coarse_flux_y(:, i, previous_j, k)) - &
            dt / dz * (coarse_flux_z(:, i, j, k) - &
              coarse_flux_z(:, i, j, previous_k))
        end do
      end do
    end do
    allocate(recovered_temperature, mold=new_hierarchy%coarse_temperature)
    call recover_reactive_temperatures_3d( &
      species, new_hierarchy%coarse_state, hierarchy%coarse_temperature, &
      distribution%coarse_count, patch%coarse_ny, patch%coarse_nz, &
      recovered_temperature, local_ok)
    if (local_ok) new_hierarchy%coarse_temperature = recovered_temperature
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    deallocate(recovered_temperature)

    call build_mpi_amr_sparse_periodic_halos_3d( &
      distribution, hierarchy%coarse_state, hierarchy%coarse_temperature, &
      coarse_start_halo, coarse_start_halo_temperature, local_ok)
    if (.not. local_ok) return
    call build_mpi_amr_sparse_periodic_halos_3d( &
      distribution, new_hierarchy%coarse_state, &
      new_hierarchy%coarse_temperature, coarse_end_halo, &
      coarse_end_halo_temperature, local_ok)
    if (.not. local_ok) return

    ratio = patch%refinement_ratio
    ratio_real = real(ratio, dp)
    subcycles = ratio * ratio
    fine_dt = dt / real(subcycles, dp)
    covered_nx = patch%coarse_i_upper - patch%coarse_i_lower + 1
    covered_ny = patch%coarse_j_upper - patch%coarse_j_lower + 1
    covered_nz = patch%coarse_k_upper - patch%coarse_k_lower + 1
    allocate(fine_x_lower(hierarchy%nvar, covered_ny, covered_nz))
    allocate(fine_x_upper(hierarchy%nvar, covered_ny, covered_nz))
    allocate(fine_y_lower(hierarchy%nvar, covered_nx, covered_nz))
    allocate(fine_y_upper(hierarchy%nvar, covered_nx, covered_nz))
    allocate(fine_z_lower(hierarchy%nvar, covered_nx, covered_ny))
    allocate(fine_z_upper(hierarchy%nvar, covered_nx, covered_ny))
    allocate(boundary_x_lower( &
      hierarchy%nvar, patch%fine_ny(), patch%fine_nz()))
    allocate(boundary_x_upper( &
      hierarchy%nvar, patch%fine_ny(), patch%fine_nz()))
    allocate(boundary_y_lower( &
      hierarchy%nvar, distribution%fine_count, patch%fine_nz()))
    allocate(boundary_y_upper( &
      hierarchy%nvar, distribution%fine_count, patch%fine_nz()))
    allocate(boundary_z_lower( &
      hierarchy%nvar, distribution%fine_count, patch%fine_ny()))
    allocate(boundary_z_upper( &
      hierarchy%nvar, distribution%fine_count, patch%fine_ny()))
    fine_x_lower = 0.0_dp
    fine_x_upper = 0.0_dp
    fine_y_lower = 0.0_dp
    fine_y_upper = 0.0_dp
    fine_z_lower = 0.0_dp
    fine_z_upper = 0.0_dp
    boundary_x_lower = 0.0_dp
    boundary_x_upper = 0.0_dp
    boundary_y_lower = 0.0_dp
    boundary_y_upper = 0.0_dp
    boundary_z_lower = 0.0_dp
    boundary_z_upper = 0.0_dp
    minimum_theta = coarse_theta

    do substep = 1, subcycles
      alpha = real(substep - 1, dp) / real(subcycles, dp)
      call compute_sparse_fine_transport_fluxes_3d( &
        distribution, species, transport, patch, coarse_start_halo, &
        coarse_start_halo_temperature, coarse_end_halo, &
        coarse_end_halo_temperature, new_hierarchy%fine_state, &
        new_hierarchy%fine_temperature, alpha, limiter, dx / ratio_real, &
        dy / ratio_real, dz / ratio_real, fine_dt, viscosity_enabled, &
        thermal_conduction_enabled, species_diffusion_enabled, &
        barodiffusion_enabled, fine_flux_x, fine_flux_y, fine_flux_z, &
        fine_theta, local_ok)
      if (.not. local_ok) return
      minimum_theta = min(minimum_theta, fine_theta)
      if (distribution%fine_count > 0) then
        do k = 1, patch%fine_nz()
          do j = 1, patch%fine_ny()
            do i = 1, distribution%fine_count
              new_hierarchy%fine_state(:, i, j, k) = &
                new_hierarchy%fine_state(:, i, j, k) - &
                fine_dt / (dx / ratio_real) * &
                  (fine_flux_x(:, i, j, k) - &
                   fine_flux_x(:, i - 1, j, k)) - &
                fine_dt / (dy / ratio_real) * &
                  (fine_flux_y(:, i, j, k) - &
                   fine_flux_y(:, i, j - 1, k)) - &
                fine_dt / (dz / ratio_real) * &
                  (fine_flux_z(:, i, j, k) - &
                   fine_flux_z(:, i, j, k - 1))
            end do
          end do
        end do
        allocate(recovered_temperature, mold=new_hierarchy%fine_temperature)
        call recover_reactive_temperatures_3d( &
          species, new_hierarchy%fine_state, &
          new_hierarchy%fine_temperature, distribution%fine_count, &
          patch%fine_ny(), patch%fine_nz(), recovered_temperature, local_ok)
        if (local_ok) &
          new_hierarchy%fine_temperature = recovered_temperature
        deallocate(recovered_temperature)
      else
        local_ok = .true.
      end if
      call sparse_collective_logical_and_3d( &
        distribution%comm, local_ok, global_ok)
      if (.not. global_ok) return
      call accumulate_sparse_fine_interface_fluxes_3d( &
        distribution, patch, fine_flux_x, fine_flux_y, fine_flux_z, &
        fine_x_lower, fine_x_upper, fine_y_lower, fine_y_upper, &
        fine_z_lower, fine_z_upper, local_ok)
      if (.not. local_ok) return
      if (distribution%fine_count > 0) then
        if (distribution%fine_first == 1) &
          boundary_x_lower = boundary_x_lower + fine_flux_x(:, 0, :, :)
        if (distribution%fine_last == patch%fine_nx()) &
          boundary_x_upper = boundary_x_upper + &
            fine_flux_x(:, distribution%fine_count, :, :)
        boundary_y_lower = boundary_y_lower + fine_flux_y(:, :, 0, :)
        boundary_y_upper = boundary_y_upper + &
          fine_flux_y(:, :, patch%fine_ny(), :)
        boundary_z_lower = boundary_z_lower + fine_flux_z(:, :, :, 0)
        boundary_z_upper = boundary_z_upper + &
          fine_flux_z(:, :, :, patch%fine_nz())
      end if
    end do
    fine_x_lower = fine_x_lower / real(subcycles, dp)
    fine_x_upper = fine_x_upper / real(subcycles, dp)
    fine_y_lower = fine_y_lower / real(subcycles, dp)
    fine_y_upper = fine_y_upper / real(subcycles, dp)
    fine_z_lower = fine_z_lower / real(subcycles, dp)
    fine_z_upper = fine_z_upper / real(subcycles, dp)
    boundary_x_lower = boundary_x_lower / real(subcycles, dp)
    boundary_x_upper = boundary_x_upper / real(subcycles, dp)
    boundary_y_lower = boundary_y_lower / real(subcycles, dp)
    boundary_y_upper = boundary_y_upper / real(subcycles, dp)
    boundary_z_lower = boundary_z_lower / real(subcycles, dp)
    boundary_z_upper = boundary_z_upper / real(subcycles, dp)

    call limit_sparse_fine_transport_interfaces_3d( &
      distribution, size(species), patch, new_hierarchy%coarse_state, &
      coarse_flux_x, coarse_flux_y, coarse_flux_z, fine_x_lower, &
      fine_x_upper, fine_y_lower, fine_y_upper, fine_z_lower, fine_z_upper, &
      boundary_x_lower, boundary_x_upper, boundary_y_lower, &
      boundary_y_upper, boundary_z_lower, boundary_z_upper, &
      new_hierarchy%fine_state, dx, dy, dz, dt, interface_theta, local_ok)
    if (.not. local_ok) return
    minimum_theta = min(minimum_theta, interface_theta)
    if (distribution%fine_count > 0) then
      allocate(recovered_temperature, mold=new_hierarchy%fine_temperature)
      call recover_reactive_temperatures_3d( &
        species, new_hierarchy%fine_state, &
        new_hierarchy%fine_temperature, distribution%fine_count, &
        patch%fine_ny(), patch%fine_nz(), recovered_temperature, local_ok)
      if (local_ok) new_hierarchy%fine_temperature = recovered_temperature
      deallocate(recovered_temperature)
    else
      local_ok = .true.
    end if
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return

    call reflux_sparse_coarse_3d( &
      distribution, patch, new_hierarchy%coarse_state, coarse_flux_x, &
      coarse_flux_y, coarse_flux_z, fine_x_lower, fine_x_upper, &
      fine_y_lower, fine_y_upper, fine_z_lower, fine_z_upper, &
      dx, dy, dz, dt, local_reflux, local_ok)
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    call MPI_Allreduce( &
      local_reflux, global_reflux, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
      distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, global_ok)
    if (.not. global_ok) return
    call average_down_sparse_hierarchy_3d( &
      distribution, patch, new_hierarchy, local_ok)
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    allocate(recovered_temperature, mold=new_hierarchy%coarse_temperature)
    call recover_reactive_temperatures_3d( &
      species, new_hierarchy%coarse_state, &
      new_hierarchy%coarse_temperature, distribution%coarse_count, &
      patch%coarse_ny, patch%coarse_nz, recovered_temperature, local_ok)
    if (local_ok) new_hierarchy%coarse_temperature = recovered_temperature
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    local_ok = new_hierarchy%is_valid(distribution, patch)
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return

    maximum_reflux_correction = global_reflux
    ok = .true.
  end subroutine advance_mpi_amr_sparse_transport_euler_3d

  subroutine limit_sparse_fine_transport_interfaces_3d( &
      distribution, nspecies, patch, coarse_state, coarse_flux_x, &
      coarse_flux_y, coarse_flux_z, fine_x_lower, fine_x_upper, &
      fine_y_lower, fine_y_upper, fine_z_lower, fine_z_upper, &
      boundary_x_lower, boundary_x_upper, boundary_y_lower, &
      boundary_y_upper, boundary_z_lower, boundary_z_upper, fine_state, &
      dx, dy, dz, dt, minimum_theta, ok)
    type(mpi_amr_sparse_distribution_3d), intent(in) :: distribution
    integer, intent(in) :: nspecies
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: coarse_state(:, :, :, :)
    real(dp), intent(in) :: coarse_flux_x(:, 0:, :, :)
    real(dp), intent(in) :: coarse_flux_y(:, :, :, :)
    real(dp), intent(in) :: coarse_flux_z(:, :, :, :)
    real(dp), intent(inout) :: fine_x_lower(:, :, :)
    real(dp), intent(inout) :: fine_x_upper(:, :, :)
    real(dp), intent(inout) :: fine_y_lower(:, :, :)
    real(dp), intent(inout) :: fine_y_upper(:, :, :)
    real(dp), intent(inout) :: fine_z_lower(:, :, :)
    real(dp), intent(inout) :: fine_z_upper(:, :, :)
    real(dp), intent(in) :: boundary_x_lower(:, :, :)
    real(dp), intent(in) :: boundary_x_upper(:, :, :)
    real(dp), intent(in) :: boundary_y_lower(:, :, :)
    real(dp), intent(in) :: boundary_y_upper(:, :, :)
    real(dp), intent(in) :: boundary_z_lower(:, :, :)
    real(dp), intent(in) :: boundary_z_upper(:, :, :)
    real(dp), intent(inout) :: fine_state(:, :, :, :)
    real(dp), intent(in) :: dx, dy, dz, dt
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok

    real(dp), allocatable :: base_state(:)
    real(dp), allocatable :: local_x_lower(:, :), local_x_upper(:, :)
    real(dp), allocatable :: local_y_lower(:, :), local_y_upper(:, :)
    real(dp), allocatable :: local_z_lower(:, :), local_z_upper(:, :)
    real(dp), allocatable :: theta_x_lower(:, :), theta_x_upper(:, :)
    real(dp), allocatable :: theta_y_lower(:, :), theta_y_upper(:, :)
    real(dp), allocatable :: theta_z_lower(:, :), theta_z_upper(:, :)
    real(dp) :: theta, fine_dx, fine_dy, fine_dz
    integer :: covered_i, global_i, local_i, local_fine_i
    integer :: fine_i, fine_j, fine_k, j, k, face_i, ratio, ierr
    integer :: covered_nx, covered_ny, covered_nz
    logical :: call_ok, global_ok, local_ok

    minimum_theta = 1.0_dp
    ok = .false.
    covered_nx = patch%coarse_i_upper - patch%coarse_i_lower + 1
    covered_ny = patch%coarse_j_upper - patch%coarse_j_lower + 1
    covered_nz = patch%coarse_k_upper - patch%coarse_k_lower + 1
    local_ok = all(ieee_is_finite([dx, dy, dz, dt]))
    if (local_ok) local_ok = dx > 0.0_dp .and. dy > 0.0_dp .and. &
      dz > 0.0_dp .and. dt > 0.0_dp
    if (local_ok) local_ok = nspecies >= 1
    if (local_ok) local_ok = patch%is_strictly_interior()
    if (local_ok) local_ok = distribution%is_valid(patch)
    if (local_ok) local_ok = all(ieee_is_finite(coarse_state))
    if (local_ok) local_ok = all(ieee_is_finite(fine_state))
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return

    ratio = patch%refinement_ratio
    fine_dx = dx / real(ratio, dp)
    fine_dy = dy / real(ratio, dp)
    fine_dz = dz / real(ratio, dp)
    allocate(base_state(size(coarse_state, 1)))
    allocate(local_x_lower(covered_ny, covered_nz), &
      local_x_upper(covered_ny, covered_nz))
    allocate(local_y_lower(covered_nx, covered_nz), &
      local_y_upper(covered_nx, covered_nz))
    allocate(local_z_lower(covered_nx, covered_ny), &
      local_z_upper(covered_nx, covered_ny))
    allocate(theta_x_lower, mold=local_x_lower)
    allocate(theta_x_upper, mold=local_x_upper)
    allocate(theta_y_lower, mold=local_y_lower)
    allocate(theta_y_upper, mold=local_y_upper)
    allocate(theta_z_lower, mold=local_z_lower)
    allocate(theta_z_upper, mold=local_z_upper)
    local_x_lower = 1.0_dp
    local_x_upper = 1.0_dp
    local_y_lower = 1.0_dp
    local_y_upper = 1.0_dp
    local_z_lower = 1.0_dp
    local_z_upper = 1.0_dp
    local_ok = .true.

    global_i = patch%coarse_i_lower - 1
    if (distribution%owns_coarse(global_i)) then
      local_i = distribution%coarse_local_index(global_i)
      face_i = global_i - distribution%coarse_first + 1
      do k = 1, covered_nz
        do j = 1, covered_ny
          base_state = coarse_state(:, local_i, &
            patch%coarse_j_lower + j - 1, patch%coarse_k_lower + k - 1) + &
            dt / dx * coarse_flux_x(:, face_i, &
              patch%coarse_j_lower + j - 1, &
              patch%coarse_k_lower + k - 1)
          call reactive_transport_interface_theta_3d( &
            nspecies, base_state, fine_x_lower(:, j, k), -dt / dx, &
            theta, call_ok)
          if (call_ok) local_x_lower(j, k) = theta
          local_ok = local_ok .and. call_ok
        end do
      end do
    end if
    global_i = patch%coarse_i_upper + 1
    if (distribution%owns_coarse(global_i)) then
      local_i = distribution%coarse_local_index(global_i)
      face_i = patch%coarse_i_upper - distribution%coarse_first + 1
      do k = 1, covered_nz
        do j = 1, covered_ny
          base_state = coarse_state(:, local_i, &
            patch%coarse_j_lower + j - 1, patch%coarse_k_lower + k - 1) - &
            dt / dx * coarse_flux_x(:, face_i, &
              patch%coarse_j_lower + j - 1, &
              patch%coarse_k_lower + k - 1)
          call reactive_transport_interface_theta_3d( &
            nspecies, base_state, fine_x_upper(:, j, k), dt / dx, &
            theta, call_ok)
          if (call_ok) local_x_upper(j, k) = theta
          local_ok = local_ok .and. call_ok
        end do
      end do
    end if

    do global_i = max(distribution%coarse_first, patch%coarse_i_lower), &
        min(distribution%coarse_last, patch%coarse_i_upper)
      local_i = distribution%coarse_local_index(global_i)
      covered_i = global_i - patch%coarse_i_lower + 1
      do k = 1, covered_nz
        base_state = coarse_state(:, local_i, patch%coarse_j_lower - 1, &
          patch%coarse_k_lower + k - 1) + &
          dt / dy * coarse_flux_y(:, local_i, &
            patch%coarse_j_lower - 1, patch%coarse_k_lower + k - 1)
        call reactive_transport_interface_theta_3d( &
          nspecies, base_state, fine_y_lower(:, covered_i, k), -dt / dy, &
          theta, call_ok)
        if (call_ok) local_y_lower(covered_i, k) = theta
        local_ok = local_ok .and. call_ok
        base_state = coarse_state(:, local_i, patch%coarse_j_upper + 1, &
          patch%coarse_k_lower + k - 1) - &
          dt / dy * coarse_flux_y(:, local_i, &
            patch%coarse_j_upper, patch%coarse_k_lower + k - 1)
        call reactive_transport_interface_theta_3d( &
          nspecies, base_state, fine_y_upper(:, covered_i, k), dt / dy, &
          theta, call_ok)
        if (call_ok) local_y_upper(covered_i, k) = theta
        local_ok = local_ok .and. call_ok
      end do
      do j = 1, covered_ny
        base_state = coarse_state(:, local_i, &
          patch%coarse_j_lower + j - 1, patch%coarse_k_lower - 1) + &
          dt / dz * coarse_flux_z(:, local_i, &
            patch%coarse_j_lower + j - 1, patch%coarse_k_lower - 1)
        call reactive_transport_interface_theta_3d( &
          nspecies, base_state, fine_z_lower(:, covered_i, j), -dt / dz, &
          theta, call_ok)
        if (call_ok) local_z_lower(covered_i, j) = theta
        local_ok = local_ok .and. call_ok
        base_state = coarse_state(:, local_i, &
          patch%coarse_j_lower + j - 1, patch%coarse_k_upper + 1) - &
          dt / dz * coarse_flux_z(:, local_i, &
            patch%coarse_j_lower + j - 1, patch%coarse_k_upper)
        call reactive_transport_interface_theta_3d( &
          nspecies, base_state, fine_z_upper(:, covered_i, j), dt / dz, &
          theta, call_ok)
        if (call_ok) local_z_upper(covered_i, j) = theta
        local_ok = local_ok .and. call_ok
      end do
    end do
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return

    call MPI_Allreduce(local_x_lower, theta_x_lower, size(local_x_lower), &
      MPI_DOUBLE_PRECISION, MPI_MIN, distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, global_ok)
    if (.not. global_ok) return
    call MPI_Allreduce(local_x_upper, theta_x_upper, size(local_x_upper), &
      MPI_DOUBLE_PRECISION, MPI_MIN, distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, global_ok)
    if (.not. global_ok) return
    call MPI_Allreduce(local_y_lower, theta_y_lower, size(local_y_lower), &
      MPI_DOUBLE_PRECISION, MPI_MIN, distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, global_ok)
    if (.not. global_ok) return
    call MPI_Allreduce(local_y_upper, theta_y_upper, size(local_y_upper), &
      MPI_DOUBLE_PRECISION, MPI_MIN, distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, global_ok)
    if (.not. global_ok) return
    call MPI_Allreduce(local_z_lower, theta_z_lower, size(local_z_lower), &
      MPI_DOUBLE_PRECISION, MPI_MIN, distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, global_ok)
    if (.not. global_ok) return
    call MPI_Allreduce(local_z_upper, theta_z_upper, size(local_z_upper), &
      MPI_DOUBLE_PRECISION, MPI_MIN, distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, global_ok)
    if (.not. global_ok) return

    local_ok = all(ieee_is_finite(theta_x_lower))
    if (local_ok) local_ok = all(ieee_is_finite(theta_x_upper))
    if (local_ok) local_ok = all(ieee_is_finite(theta_y_lower))
    if (local_ok) local_ok = all(ieee_is_finite(theta_y_upper))
    if (local_ok) local_ok = all(ieee_is_finite(theta_z_lower))
    if (local_ok) local_ok = all(ieee_is_finite(theta_z_upper))
    if (local_ok) local_ok = minval(theta_x_lower) >= 0.0_dp
    if (local_ok) local_ok = maxval(theta_x_lower) <= 1.0_dp
    if (local_ok) local_ok = minval(theta_x_upper) >= 0.0_dp
    if (local_ok) local_ok = maxval(theta_x_upper) <= 1.0_dp
    if (local_ok) local_ok = minval(theta_y_lower) >= 0.0_dp
    if (local_ok) local_ok = maxval(theta_y_lower) <= 1.0_dp
    if (local_ok) local_ok = minval(theta_y_upper) >= 0.0_dp
    if (local_ok) local_ok = maxval(theta_y_upper) <= 1.0_dp
    if (local_ok) local_ok = minval(theta_z_lower) >= 0.0_dp
    if (local_ok) local_ok = maxval(theta_z_lower) <= 1.0_dp
    if (local_ok) local_ok = minval(theta_z_upper) >= 0.0_dp
    if (local_ok) local_ok = maxval(theta_z_upper) <= 1.0_dp
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return

    do k = 1, covered_nz
      do j = 1, covered_ny
        theta = theta_x_lower(j, k)
        if (theta < 1.0_dp) then
          fine_x_lower(:, j, k) = theta * fine_x_lower(:, j, k)
          if (distribution%fine_count > 0 .and. &
              distribution%fine_first == 1) then
            do fine_k = (k - 1) * ratio + 1, k * ratio
              do fine_j = (j - 1) * ratio + 1, j * ratio
                fine_state(:, 1, fine_j, fine_k) = &
                  fine_state(:, 1, fine_j, fine_k) + &
                  dt / fine_dx * (theta - 1.0_dp) * &
                    boundary_x_lower(:, fine_j, fine_k)
              end do
            end do
          end if
        end if
        theta = theta_x_upper(j, k)
        if (theta < 1.0_dp) then
          fine_x_upper(:, j, k) = theta * fine_x_upper(:, j, k)
          if (distribution%fine_count > 0 .and. &
              distribution%fine_last == patch%fine_nx()) then
            do fine_k = (k - 1) * ratio + 1, k * ratio
              do fine_j = (j - 1) * ratio + 1, j * ratio
                fine_state(:, distribution%fine_count, fine_j, fine_k) = &
                  fine_state(:, distribution%fine_count, fine_j, fine_k) - &
                  dt / fine_dx * (theta - 1.0_dp) * &
                    boundary_x_upper(:, fine_j, fine_k)
              end do
            end do
          end if
        end if
      end do
    end do

    do k = 1, covered_nz
      do covered_i = 1, covered_nx
        theta = theta_y_lower(covered_i, k)
        if (theta < 1.0_dp) &
          fine_y_lower(:, covered_i, k) = &
            theta * fine_y_lower(:, covered_i, k)
        theta = theta_y_upper(covered_i, k)
        if (theta < 1.0_dp) &
          fine_y_upper(:, covered_i, k) = &
            theta * fine_y_upper(:, covered_i, k)
      end do
    end do
    do local_fine_i = 1, distribution%fine_count
      fine_i = distribution%fine_first + local_fine_i - 1
      covered_i = 1 + (fine_i - 1) / ratio
      do fine_k = 1, patch%fine_nz()
        k = 1 + (fine_k - 1) / ratio
        theta = theta_y_lower(covered_i, k)
        if (theta < 1.0_dp) &
          fine_state(:, local_fine_i, 1, fine_k) = &
            fine_state(:, local_fine_i, 1, fine_k) + &
            dt / fine_dy * (theta - 1.0_dp) * &
              boundary_y_lower(:, local_fine_i, fine_k)
        theta = theta_y_upper(covered_i, k)
        if (theta < 1.0_dp) &
          fine_state(:, local_fine_i, patch%fine_ny(), fine_k) = &
            fine_state(:, local_fine_i, patch%fine_ny(), fine_k) - &
            dt / fine_dy * (theta - 1.0_dp) * &
              boundary_y_upper(:, local_fine_i, fine_k)
      end do
    end do

    do j = 1, covered_ny
      do covered_i = 1, covered_nx
        theta = theta_z_lower(covered_i, j)
        if (theta < 1.0_dp) &
          fine_z_lower(:, covered_i, j) = &
            theta * fine_z_lower(:, covered_i, j)
        theta = theta_z_upper(covered_i, j)
        if (theta < 1.0_dp) &
          fine_z_upper(:, covered_i, j) = &
            theta * fine_z_upper(:, covered_i, j)
      end do
    end do
    do local_fine_i = 1, distribution%fine_count
      fine_i = distribution%fine_first + local_fine_i - 1
      covered_i = 1 + (fine_i - 1) / ratio
      do fine_j = 1, patch%fine_ny()
        j = 1 + (fine_j - 1) / ratio
        theta = theta_z_lower(covered_i, j)
        if (theta < 1.0_dp) &
          fine_state(:, local_fine_i, fine_j, 1) = &
            fine_state(:, local_fine_i, fine_j, 1) + &
            dt / fine_dz * (theta - 1.0_dp) * &
              boundary_z_lower(:, local_fine_i, fine_j)
        theta = theta_z_upper(covered_i, j)
        if (theta < 1.0_dp) &
          fine_state(:, local_fine_i, fine_j, patch%fine_nz()) = &
            fine_state(:, local_fine_i, fine_j, patch%fine_nz()) - &
            dt / fine_dz * (theta - 1.0_dp) * &
              boundary_z_upper(:, local_fine_i, fine_j)
      end do
    end do

    local_ok = all(ieee_is_finite(fine_x_lower))
    if (local_ok) local_ok = all(ieee_is_finite(fine_x_upper))
    if (local_ok) local_ok = all(ieee_is_finite(fine_y_lower))
    if (local_ok) local_ok = all(ieee_is_finite(fine_y_upper))
    if (local_ok) local_ok = all(ieee_is_finite(fine_z_lower))
    if (local_ok) local_ok = all(ieee_is_finite(fine_z_upper))
    if (local_ok) local_ok = all(ieee_is_finite(fine_state))
    if (local_ok) then
      minimum_theta = min( &
        minval(theta_x_lower), minval(theta_x_upper), &
        minval(theta_y_lower), minval(theta_y_upper), &
        minval(theta_z_lower), minval(theta_z_upper))
      local_ok = ieee_is_finite(minimum_theta)
      if (local_ok) local_ok = minimum_theta >= 0.0_dp
      if (local_ok) local_ok = minimum_theta <= 1.0_dp
    end if
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, ok)
  end subroutine limit_sparse_fine_transport_interfaces_3d

  subroutine compute_sparse_periodic_transport_fluxes_3d( &
      distribution, species, transport, patch, state, temperature, &
      dx, dy, dz, dt, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, face_flux_x, &
      face_flux_y, face_flux_z, minimum_theta, ok)
    type(mpi_amr_sparse_distribution_3d), intent(in) :: distribution
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    real(dp), intent(in) :: dx, dy, dz, dt
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), allocatable, intent(out) :: face_flux_x(:, :, :, :)
    real(dp), allocatable, intent(out) :: face_flux_y(:, :, :, :)
    real(dp), allocatable, intent(out) :: face_flux_z(:, :, :, :)
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok

    real(dp), allocatable :: halo_state(:, :, :, :)
    real(dp), allocatable :: halo_temperature(:, :, :)
    real(dp), allocatable :: ghost_state(:, :, :, :)
    real(dp), allocatable :: ghost_temperature(:, :, :)
    real(dp), allocatable :: explicit_flux_x(:, :, :, :)
    real(dp), allocatable :: explicit_flux_y(:, :, :, :)
    real(dp), allocatable :: explicit_flux_z(:, :, :, :)
    real(dp), allocatable :: theta_cell(:, :, :)
    real(dp), allocatable :: theta_lower(:, :), theta_upper(:, :)
    real(dp) :: local_theta
    logical :: local_ok, global_ok
    integer :: i, j, k, halo_i, source_j, source_k
    integer :: local_count, ny, nz, nvar, left_rank, right_rank, ierr

    minimum_theta = 1.0_dp
    ok = .false.
    local_count = distribution%coarse_count
    ny = patch%coarse_ny
    nz = patch%coarse_nz
    nvar = reactive_nvar(size(species))
    call build_mpi_amr_sparse_periodic_halos_3d( &
      distribution, state, temperature, halo_state, halo_temperature, local_ok)
    if (.not. local_ok) return

    allocate(ghost_state(nvar, 0:local_count + 1, 0:ny + 1, 0:nz + 1))
    allocate(ghost_temperature(0:local_count + 1, 0:ny + 1, 0:nz + 1))
    do k = 0, nz + 1
      source_k = 1 + modulo(k - 1, nz)
      do j = 0, ny + 1
        source_j = 1 + modulo(j - 1, ny)
        do i = 0, local_count + 1
          halo_i = i + mpi_amr_sparse_ghost_width_3d
          ghost_state(:, i, j, k) = &
            halo_state(:, halo_i, source_j, source_k)
          ghost_temperature(i, j, k) = &
            halo_temperature(halo_i, source_j, source_k)
        end do
      end do
    end do
    allocate(explicit_flux_x(nvar, 0:local_count, ny, nz))
    allocate(explicit_flux_y(nvar, local_count, 0:ny, nz))
    allocate(explicit_flux_z(nvar, local_count, ny, 0:nz))
    allocate(theta_cell(local_count, ny, nz))
    call reactive_transport_ghosted_fluxes_3d( &
      species, transport, ghost_state, ghost_temperature, &
      local_count, ny, nz, dx, dy, dz, dt, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, explicit_flux_x, explicit_flux_y, &
      explicit_flux_z, local_theta, local_ok, &
      theta_cell_output=theta_cell)
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return

    allocate(theta_lower(ny, nz), theta_upper(ny, nz))
    left_rank = modulo(distribution%rank - 1, distribution%nranks)
    right_rank = modulo(distribution%rank + 1, distribution%nranks)
    call exchange_sparse_theta_x_3d( &
      distribution%comm, left_rank, right_rank, theta_cell, &
      theta_lower, theta_upper, local_ok)
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    call reactive_transport_ghosted_fluxes_3d( &
      species, transport, ghost_state, ghost_temperature, &
      local_count, ny, nz, dx, dy, dz, dt, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, explicit_flux_x, explicit_flux_y, &
      explicit_flux_z, local_theta, local_ok, &
      exterior_theta_x_lower=theta_lower, &
      exterior_theta_x_upper=theta_upper, periodic_theta_y=.true., &
      periodic_theta_z=.true.)
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return

    allocate(face_flux_x(nvar, 0:local_count, ny, nz))
    allocate(face_flux_y(nvar, local_count, ny, nz))
    allocate(face_flux_z(nvar, local_count, ny, nz))
    face_flux_x = explicit_flux_x
    face_flux_y = explicit_flux_y(:, :, 1:ny, :)
    face_flux_z = explicit_flux_z(:, :, :, 1:nz)
    call MPI_Allreduce( &
      local_theta, minimum_theta, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
      distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, global_ok)
    if (.not. global_ok) return
    ok = ieee_is_finite(minimum_theta)
    if (ok) ok = minimum_theta >= 0.0_dp
    if (ok) ok = minimum_theta <= 1.0_dp
  end subroutine compute_sparse_periodic_transport_fluxes_3d

  subroutine compute_sparse_fine_transport_fluxes_3d( &
      distribution, species, transport, patch, coarse_start, &
      coarse_start_temperature, coarse_end, coarse_end_temperature, &
      state, temperature, alpha, limiter, dx, dy, dz, dt, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, barodiffusion_enabled, face_flux_x, &
      face_flux_y, face_flux_z, minimum_theta, ok)
    type(mpi_amr_sparse_distribution_3d), intent(in) :: distribution
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: coarse_start(:, :, :, :)
    real(dp), intent(in) :: coarse_start_temperature(:, :, :)
    real(dp), intent(in) :: coarse_end(:, :, :, :)
    real(dp), intent(in) :: coarse_end_temperature(:, :, :)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    real(dp), intent(in) :: alpha, dx, dy, dz, dt
    character(len=*), intent(in) :: limiter
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled, barodiffusion_enabled
    real(dp), allocatable, intent(out) :: face_flux_x(:, :, :, :)
    real(dp), allocatable, intent(out) :: face_flux_y(:, :, :, :)
    real(dp), allocatable, intent(out) :: face_flux_z(:, :, :, :)
    real(dp), intent(out) :: minimum_theta
    logical, intent(out) :: ok

    real(dp), allocatable :: extended_state(:, :, :, :)
    real(dp), allocatable :: extended_temperature(:, :, :)
    real(dp), allocatable :: ghost_state(:, :, :, :)
    real(dp), allocatable :: ghost_temperature(:, :, :)
    real(dp), allocatable :: theta_cell(:, :, :)
    real(dp), allocatable :: theta_lower(:, :), theta_upper(:, :)
    real(dp) :: local_theta
    logical :: local_ok, global_ok
    integer :: local_count, ny, nz, nvar, ierr

    minimum_theta = 1.0_dp
    ok = .false.
    local_count = distribution%fine_count
    ny = patch%fine_ny()
    nz = patch%fine_nz()
    nvar = reactive_nvar(size(species))
    local_ok = all(ieee_is_finite([alpha, dx, dy, dz, dt]))
    if (local_ok) local_ok = alpha >= 0.0_dp .and. alpha <= 1.0_dp
    if (local_ok) local_ok = dx > 0.0_dp .and. dy > 0.0_dp .and. &
      dz > 0.0_dp .and. dt > 0.0_dp
    if (local_ok) local_ok = size(state, 1) == nvar .and. &
      size(state, 2) == local_count .and. size(state, 3) == ny .and. &
      size(state, 4) == nz
    if (local_ok) local_ok = all(shape(temperature) == [local_count, ny, nz])
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return

    allocate(face_flux_x(nvar, 0:local_count, ny, nz))
    allocate(face_flux_y(nvar, local_count, 0:ny, nz))
    allocate(face_flux_z(nvar, local_count, ny, 0:nz))
    face_flux_x = 0.0_dp
    face_flux_y = 0.0_dp
    face_flux_z = 0.0_dp
    if (local_count > 0) then
      call build_sparse_fine_extended_state_3d( &
        distribution, species, patch, coarse_start, &
        coarse_start_temperature, coarse_end, coarse_end_temperature, &
        state, temperature, alpha, "characteristic_plm", limiter, &
        extended_state, extended_temperature, local_ok)
    else
      local_ok = .true.
    end if
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return

    if (local_count > 0) then
      allocate(ghost_state(nvar, 0:local_count + 1, 0:ny + 1, 0:nz + 1))
      allocate(ghost_temperature(0:local_count + 1, 0:ny + 1, 0:nz + 1))
      ghost_state = extended_state(:, 2:local_count + 3, 2:ny + 3, 2:nz + 3)
      ghost_temperature = &
        extended_temperature(2:local_count + 3, 2:ny + 3, 2:nz + 3)
      allocate(theta_cell(local_count, ny, nz))
      call reactive_transport_ghosted_fluxes_3d( &
        species, transport, ghost_state, ghost_temperature, &
        local_count, ny, nz, dx, dy, dz, dt, viscosity_enabled, &
        thermal_conduction_enabled, species_diffusion_enabled, &
        barodiffusion_enabled, face_flux_x, face_flux_y, face_flux_z, &
        local_theta, local_ok, theta_cell_output=theta_cell)
    else
      local_theta = 1.0_dp
      local_ok = .true.
    end if
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return

    if (local_count > 0) then
      allocate(theta_lower(ny, nz), theta_upper(ny, nz))
      call exchange_sparse_theta_x_3d( &
        distribution%comm, distribution%previous_fine_rank, &
        distribution%next_fine_rank, theta_cell, theta_lower, &
        theta_upper, local_ok)
    else
      local_ok = .true.
    end if
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return

    if (local_count > 0) then
      if (distribution%previous_fine_rank /= MPI_PROC_NULL .and. &
          distribution%next_fine_rank /= MPI_PROC_NULL) then
        call reactive_transport_ghosted_fluxes_3d( &
          species, transport, ghost_state, ghost_temperature, &
          local_count, ny, nz, dx, dy, dz, dt, viscosity_enabled, &
          thermal_conduction_enabled, species_diffusion_enabled, &
          barodiffusion_enabled, face_flux_x, face_flux_y, face_flux_z, &
          local_theta, local_ok, exterior_theta_x_lower=theta_lower, &
          exterior_theta_x_upper=theta_upper)
      else if (distribution%previous_fine_rank /= MPI_PROC_NULL) then
        call reactive_transport_ghosted_fluxes_3d( &
          species, transport, ghost_state, ghost_temperature, &
          local_count, ny, nz, dx, dy, dz, dt, viscosity_enabled, &
          thermal_conduction_enabled, species_diffusion_enabled, &
          barodiffusion_enabled, face_flux_x, face_flux_y, face_flux_z, &
          local_theta, local_ok, exterior_theta_x_lower=theta_lower)
      else if (distribution%next_fine_rank /= MPI_PROC_NULL) then
        call reactive_transport_ghosted_fluxes_3d( &
          species, transport, ghost_state, ghost_temperature, &
          local_count, ny, nz, dx, dy, dz, dt, viscosity_enabled, &
          thermal_conduction_enabled, species_diffusion_enabled, &
          barodiffusion_enabled, face_flux_x, face_flux_y, face_flux_z, &
          local_theta, local_ok, exterior_theta_x_upper=theta_upper)
      else
        call reactive_transport_ghosted_fluxes_3d( &
          species, transport, ghost_state, ghost_temperature, &
          local_count, ny, nz, dx, dy, dz, dt, viscosity_enabled, &
          thermal_conduction_enabled, species_diffusion_enabled, &
          barodiffusion_enabled, face_flux_x, face_flux_y, face_flux_z, &
          local_theta, local_ok)
      end if
    else
      local_theta = 1.0_dp
      local_ok = .true.
    end if
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    call MPI_Allreduce( &
      local_theta, minimum_theta, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
      distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, global_ok)
    if (.not. global_ok) return
    ok = ieee_is_finite(minimum_theta)
    if (ok) ok = minimum_theta >= 0.0_dp
    if (ok) ok = minimum_theta <= 1.0_dp
  end subroutine compute_sparse_fine_transport_fluxes_3d

  subroutine exchange_sparse_theta_x_3d( &
      comm, left_rank, right_rank, theta_cell, theta_lower, theta_upper, ok)
    type(MPI_Comm), intent(in) :: comm
    integer, intent(in) :: left_rank, right_rank
    real(dp), intent(in) :: theta_cell(:, :, :)
    real(dp), intent(out) :: theta_lower(:, :), theta_upper(:, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: send_values(:), receive_values(:)
    integer :: payload_size, ierr
    type(MPI_Status) :: status

    ok = .false.
    theta_lower = 1.0_dp
    theta_upper = 1.0_dp
    if (size(theta_cell, 1) < 1) return
    if (any(shape(theta_lower) /= shape(theta_cell(1, :, :))) .or. &
        any(shape(theta_upper) /= shape(theta_cell(1, :, :)))) return
    payload_size = size(theta_lower)
    allocate(send_values(payload_size), receive_values(payload_size))
    send_values = reshape(theta_cell(1, :, :), [payload_size])
    receive_values = 1.0_dp
    call MPI_Sendrecv( &
      send_values, payload_size, MPI_DOUBLE_PRECISION, left_rank, &
      sparse_theta_left_tag_3d, receive_values, payload_size, &
      MPI_DOUBLE_PRECISION, right_rank, sparse_theta_left_tag_3d, &
      comm, status, ierr)
    if (ierr /= MPI_SUCCESS) return
    if (right_rank /= MPI_PROC_NULL) &
      theta_upper = reshape(receive_values, shape(theta_upper))

    send_values = reshape( &
      theta_cell(size(theta_cell, 1), :, :), [payload_size])
    receive_values = 1.0_dp
    call MPI_Sendrecv( &
      send_values, payload_size, MPI_DOUBLE_PRECISION, right_rank, &
      sparse_theta_right_tag_3d, receive_values, payload_size, &
      MPI_DOUBLE_PRECISION, left_rank, sparse_theta_right_tag_3d, &
      comm, status, ierr)
    if (ierr /= MPI_SUCCESS) return
    if (left_rank /= MPI_PROC_NULL) &
      theta_lower = reshape(receive_values, shape(theta_lower))
    ok = all(ieee_is_finite(theta_lower))
    if (ok) ok = all(ieee_is_finite(theta_upper))
    if (ok) ok = minval(theta_lower) >= 0.0_dp
    if (ok) ok = maxval(theta_lower) <= 1.0_dp
    if (ok) ok = minval(theta_upper) >= 0.0_dp
    if (ok) ok = maxval(theta_upper) <= 1.0_dp
  end subroutine exchange_sparse_theta_x_3d

  subroutine advance_mpi_amr_sparse_hydro_3d( &
      distribution, species, patch, hierarchy, dx, dy, dz, dt, &
      riemann_solver, reconstruction, limiter, maximum_reflux_correction, ok)
    type(mpi_amr_sparse_distribution_3d), intent(in) :: distribution
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    type(mpi_amr_sparse_hierarchy_3d), intent(inout) :: hierarchy
    real(dp), intent(in) :: dx, dy, dz, dt
    character(len=*), intent(in) :: riemann_solver, reconstruction, limiter
    real(dp), intent(out) :: maximum_reflux_correction
    logical, intent(out) :: ok

    type(mpi_amr_sparse_hierarchy_3d) :: candidate
    real(dp), allocatable :: coarse_start(:, :, :, :)
    real(dp), allocatable :: coarse_start_temperature(:, :, :)
    real(dp), allocatable :: coarse_start_halo(:, :, :, :)
    real(dp), allocatable :: coarse_start_halo_temperature(:, :, :)
    real(dp), allocatable :: coarse_end_halo(:, :, :, :)
    real(dp), allocatable :: coarse_end_halo_temperature(:, :, :)
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
    real(dp) :: alpha_start, alpha_end, fine_dt
    real(dp) :: local_reflux, global_reflux
    logical :: local_ok, global_ok
    integer :: covered_nx, covered_ny, covered_nz
    integer :: ratio, substep, ierr

    maximum_reflux_correction = 0.0_dp
    ok = .false.
    local_ok = all(ieee_is_finite([dx, dy, dz, dt]))
    if (local_ok) local_ok = dx > 0.0_dp .and. dy > 0.0_dp .and. &
      dz > 0.0_dp .and. dt > 0.0_dp
    if (local_ok) local_ok = hierarchy%nvar == reactive_nvar(size(species))
    if (local_ok) local_ok = hierarchy%is_valid(distribution, patch)
    if (local_ok) local_ok = trim(reconstruction) == "pcm" .or. &
      trim(reconstruction) == "characteristic_plm"
    if (local_ok) local_ok = trim(limiter) == "minmod" .or. &
      trim(limiter) == "mc"
    call sparse_collective_contract_matches_3d( &
      distribution, species, patch, [dx, dy, dz, dt], &
      trim(riemann_solver) // "|" // trim(reconstruction) // "|" // &
        trim(limiter), local_ok, global_ok)
    if (.not. global_ok) return

    ratio = patch%refinement_ratio
    covered_nx = patch%coarse_i_upper - patch%coarse_i_lower + 1
    covered_ny = patch%coarse_j_upper - patch%coarse_j_lower + 1
    covered_nz = patch%coarse_k_upper - patch%coarse_k_lower + 1
    candidate = hierarchy
    coarse_start = hierarchy%coarse_state
    allocate(coarse_start_temperature, mold=hierarchy%coarse_temperature)
    call recover_reactive_temperatures_3d( &
      species, coarse_start, hierarchy%coarse_temperature, &
      distribution%coarse_count, patch%coarse_ny, patch%coarse_nz, &
      coarse_start_temperature, local_ok)
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    candidate%coarse_temperature = coarse_start_temperature
    call advance_mpi_amr_sparse_periodic_hydro_3d( &
      distribution, species, patch, candidate%coarse_state, &
      candidate%coarse_temperature, dx, dy, dz, dt, riemann_solver, &
      reconstruction, limiter, coarse_flux_x, coarse_flux_y, &
      coarse_flux_z, local_ok)
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    call build_mpi_amr_sparse_periodic_halos_3d( &
      distribution, coarse_start, coarse_start_temperature, &
      coarse_start_halo, coarse_start_halo_temperature, local_ok)
    if (.not. local_ok) return
    call build_mpi_amr_sparse_periodic_halos_3d( &
      distribution, candidate%coarse_state, candidate%coarse_temperature, &
      coarse_end_halo, coarse_end_halo_temperature, local_ok)
    if (.not. local_ok) return

    allocate(fine_x_lower(hierarchy%nvar, covered_ny, covered_nz))
    allocate(fine_x_upper(hierarchy%nvar, covered_ny, covered_nz))
    allocate(fine_y_lower(hierarchy%nvar, covered_nx, covered_nz))
    allocate(fine_y_upper(hierarchy%nvar, covered_nx, covered_nz))
    allocate(fine_z_lower(hierarchy%nvar, covered_nx, covered_ny))
    allocate(fine_z_upper(hierarchy%nvar, covered_nx, covered_ny))
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
      call advance_sparse_fine_patch_ssprk2_3d( &
        distribution, species, patch, coarse_start_halo, &
        coarse_start_halo_temperature, coarse_end_halo, &
        coarse_end_halo_temperature, candidate%fine_state, &
        candidate%fine_temperature, dx / real(ratio, dp), &
        dy / real(ratio, dp), dz / real(ratio, dp), fine_dt, &
        alpha_start, alpha_end, riemann_solver, reconstruction, limiter, &
        fine_flux_x, fine_flux_y, fine_flux_z, local_ok)
      call sparse_collective_logical_and_3d( &
        distribution%comm, local_ok, global_ok)
      if (.not. global_ok) return
      call accumulate_sparse_fine_interface_fluxes_3d( &
        distribution, patch, fine_flux_x, fine_flux_y, fine_flux_z, &
        fine_x_lower, fine_x_upper, fine_y_lower, fine_y_upper, &
        fine_z_lower, fine_z_upper, local_ok)
      call sparse_collective_logical_and_3d( &
        distribution%comm, local_ok, global_ok)
      if (.not. global_ok) return
    end do
    fine_x_lower = fine_x_lower / real(ratio, dp)
    fine_x_upper = fine_x_upper / real(ratio, dp)
    fine_y_lower = fine_y_lower / real(ratio, dp)
    fine_y_upper = fine_y_upper / real(ratio, dp)
    fine_z_lower = fine_z_lower / real(ratio, dp)
    fine_z_upper = fine_z_upper / real(ratio, dp)

    call reflux_sparse_coarse_3d( &
      distribution, patch, candidate%coarse_state, coarse_flux_x, &
      coarse_flux_y, coarse_flux_z, fine_x_lower, fine_x_upper, &
      fine_y_lower, fine_y_upper, fine_z_lower, fine_z_upper, &
      dx, dy, dz, dt, local_reflux, local_ok)
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    call MPI_Allreduce( &
      local_reflux, global_reflux, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
      distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, global_ok)
    if (.not. global_ok) return
    call average_down_sparse_hierarchy_3d( &
      distribution, patch, candidate, local_ok)
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    allocate(synchronized_temperature, mold=candidate%coarse_temperature)
    call recover_reactive_temperatures_3d( &
      species, candidate%coarse_state, candidate%coarse_temperature, &
      distribution%coarse_count, patch%coarse_ny, patch%coarse_nz, &
      synchronized_temperature, local_ok)
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    candidate%coarse_temperature = synchronized_temperature
    local_ok = candidate%is_valid(distribution, patch)
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return

    hierarchy = candidate
    maximum_reflux_correction = global_reflux
    ok = .true.
  end subroutine advance_mpi_amr_sparse_hydro_3d

  subroutine advance_mpi_amr_sparse_periodic_hydro_3d( &
      distribution, species, patch, state, temperature, dx, dy, dz, dt, &
      riemann_solver, reconstruction, limiter, face_flux_x, &
      face_flux_y, face_flux_z, ok)
    type(mpi_amr_sparse_distribution_3d), intent(in) :: distribution
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(inout) :: state(:, :, :, :), temperature(:, :, :)
    real(dp), intent(in) :: dx, dy, dz, dt
    character(len=*), intent(in) :: riemann_solver, reconstruction, limiter
    real(dp), allocatable, intent(out) :: face_flux_x(:, :, :, :)
    real(dp), allocatable, intent(out) :: face_flux_y(:, :, :, :)
    real(dp), allocatable, intent(out) :: face_flux_z(:, :, :, :)
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
    integer :: nvar, local_count, ny, nz

    ok = .false.
    nvar = reactive_nvar(size(species))
    local_count = size(state, 2)
    ny = size(state, 3)
    nz = size(state, 4)
    local_ok = all(ieee_is_finite([dx, dy, dz, dt]))
    if (local_ok) local_ok = dx > 0.0_dp .and. dy > 0.0_dp .and. &
      dz > 0.0_dp .and. dt > 0.0_dp
    if (local_ok) local_ok = distribution%is_valid(patch)
    if (local_ok) local_ok = size(state, 1) == nvar .and. &
      local_count == distribution%coarse_count .and. &
      ny == patch%coarse_ny .and. nz == patch%coarse_nz
    if (local_ok) local_ok = all(shape(temperature) == [local_count, ny, nz])
    if (local_ok) local_ok = all(ieee_is_finite(state))
    if (local_ok) local_ok = all(ieee_is_finite(temperature))
    if (local_ok) local_ok = minval(temperature) > 0.0_dp
    if (local_ok) local_ok = trim(reconstruction) == "pcm" .or. &
      trim(reconstruction) == "characteristic_plm"
    if (local_ok) local_ok = trim(limiter) == "minmod" .or. &
      trim(limiter) == "mc"
    call sparse_collective_contract_matches_3d( &
      distribution, species, patch, [dx, dy, dz, dt], &
      trim(riemann_solver) // "|" // trim(reconstruction) // "|" // &
        trim(limiter), local_ok, global_ok)
    if (.not. global_ok) return

    old_state = state
    allocate(old_temperature(local_count, ny, nz))
    call recover_reactive_temperatures_3d( &
      species, old_state, temperature, local_count, ny, nz, &
      old_temperature, local_ok)
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    call compute_sparse_periodic_face_fluxes_3d( &
      distribution, species, old_state, old_temperature, &
      riemann_solver, reconstruction, limiter, first_flux_x, &
      first_flux_y, first_flux_z, local_ok)
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    call update_sparse_periodic_stage_3d( &
      old_state, old_state, first_flux_x, first_flux_y, first_flux_z, &
      dx, dy, dz, dt, .false., stage_state, local_ok)
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    allocate(stage_temperature(local_count, ny, nz))
    call recover_reactive_temperatures_3d( &
      species, stage_state, old_temperature, local_count, ny, nz, &
      stage_temperature, local_ok)
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return

    call compute_sparse_periodic_face_fluxes_3d( &
      distribution, species, stage_state, stage_temperature, &
      riemann_solver, reconstruction, limiter, second_flux_x, &
      second_flux_y, second_flux_z, local_ok)
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    call update_sparse_periodic_stage_3d( &
      old_state, stage_state, second_flux_x, second_flux_y, second_flux_z, &
      dx, dy, dz, dt, .true., candidate_state, local_ok)
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    allocate(candidate_temperature(local_count, ny, nz))
    call recover_reactive_temperatures_3d( &
      species, candidate_state, stage_temperature, local_count, ny, nz, &
      candidate_temperature, local_ok)
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return

    allocate(face_flux_x(nvar, 0:local_count, ny, nz))
    allocate(face_flux_y(nvar, local_count, ny, nz))
    allocate(face_flux_z(nvar, local_count, ny, nz))
    face_flux_x = 0.5_dp * (first_flux_x + second_flux_x)
    face_flux_y = 0.5_dp * (first_flux_y + second_flux_y)
    face_flux_z = 0.5_dp * (first_flux_z + second_flux_z)
    state = candidate_state
    temperature = candidate_temperature
    ok = .true.
  end subroutine advance_mpi_amr_sparse_periodic_hydro_3d

  subroutine compute_sparse_periodic_face_fluxes_3d( &
      distribution, species, state, temperature, riemann_solver, &
      reconstruction, limiter, face_flux_x, face_flux_y, face_flux_z, ok)
    type(mpi_amr_sparse_distribution_3d), intent(in) :: distribution
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    character(len=*), intent(in) :: riemann_solver, reconstruction, limiter
    real(dp), allocatable, intent(out) :: face_flux_x(:, :, :, :)
    real(dp), allocatable, intent(out) :: face_flux_y(:, :, :, :)
    real(dp), allocatable, intent(out) :: face_flux_z(:, :, :, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: halo_state(:, :, :, :)
    real(dp), allocatable :: halo_temperature(:, :, :)
    real(dp), allocatable :: extended_flux_x(:, :, :, :)
    real(dp), allocatable :: extended_flux_y(:, :, :, :)
    real(dp), allocatable :: extended_flux_z(:, :, :, :)
    logical :: face_ok
    integer :: i, j, k, next_j, next_k
    integer :: nvar, local_count, ny, nz

    ok = .false.
    nvar = size(state, 1)
    local_count = size(state, 2)
    ny = size(state, 3)
    nz = size(state, 4)
    allocate(face_flux_x(nvar, 0:local_count, ny, nz))
    allocate(face_flux_y(nvar, local_count, ny, nz))
    allocate(face_flux_z(nvar, local_count, ny, nz))
    face_flux_x = 0.0_dp
    face_flux_y = 0.0_dp
    face_flux_z = 0.0_dp
    call build_mpi_amr_sparse_periodic_halos_3d( &
      distribution, state, temperature, halo_state, halo_temperature, face_ok)
    if (.not. face_ok) return

    select case (trim(reconstruction))
    case ("characteristic_plm")
      allocate(extended_flux_x, mold=halo_state)
      allocate(extended_flux_y, mold=halo_state)
      allocate(extended_flux_z, mold=halo_state)
      call compute_reactive_plm_slab_face_fluxes_3d( &
        species, halo_state, halo_temperature, local_count + 4, ny, nz, &
        2, local_count + 2, limiter, riemann_solver, &
        extended_flux_x, extended_flux_y, extended_flux_z, face_ok)
      if (.not. face_ok) return
      do k = 1, nz
        do j = 1, ny
          do i = 0, local_count
            face_flux_x(:, i, j, k) = &
              extended_flux_x(:, i + 2, j, k)
          end do
          do i = 1, local_count
            face_flux_y(:, i, j, k) = &
              extended_flux_y(:, i + 2, j, k)
            face_flux_z(:, i, j, k) = &
              extended_flux_z(:, i + 2, j, k)
          end do
        end do
      end do
    case ("pcm")
      outer: do k = 1, nz
        next_k = modulo(k, nz) + 1
        do j = 1, ny
          next_j = modulo(j, ny) + 1
          do i = 0, local_count
            call reactive_riemann_flux_x( &
              species, halo_state(:, i + 2, j, k), &
              halo_state(:, i + 3, j, k), &
              halo_temperature(i + 2, j, k), &
              halo_temperature(i + 3, j, k), riemann_solver, &
              face_flux_x(:, i, j, k), face_ok)
            if (.not. face_ok) exit outer
          end do
          do i = 1, local_count
            call reactive_riemann_flux_y( &
              species, halo_state(:, i + 2, j, k), &
              halo_state(:, i + 2, next_j, k), &
              halo_temperature(i + 2, j, k), &
              halo_temperature(i + 2, next_j, k), riemann_solver, &
              face_flux_y(:, i, j, k), face_ok)
            if (.not. face_ok) exit outer
            call reactive_riemann_flux_z( &
              species, halo_state(:, i + 2, j, k), &
              halo_state(:, i + 2, j, next_k), &
              halo_temperature(i + 2, j, k), &
              halo_temperature(i + 2, j, next_k), riemann_solver, &
              face_flux_z(:, i, j, k), face_ok)
            if (.not. face_ok) exit outer
          end do
        end do
      end do outer
      if (.not. face_ok) return
    case default
      return
    end select
    ok = all(ieee_is_finite(face_flux_x)) .and. &
      all(ieee_is_finite(face_flux_y)) .and. &
      all(ieee_is_finite(face_flux_z))
  end subroutine compute_sparse_periodic_face_fluxes_3d

  subroutine update_sparse_periodic_stage_3d( &
      old_state, stage_base, face_flux_x, face_flux_y, face_flux_z, &
      dx, dy, dz, dt, second_stage, candidate, ok)
    real(dp), intent(in) :: old_state(:, :, :, :)
    real(dp), intent(in) :: stage_base(:, :, :, :)
    real(dp), intent(in) :: face_flux_x(:, 0:, :, :)
    real(dp), intent(in) :: face_flux_y(:, :, :, :)
    real(dp), intent(in) :: face_flux_z(:, :, :, :)
    real(dp), intent(in) :: dx, dy, dz, dt
    logical, intent(in) :: second_stage
    real(dp), allocatable, intent(out) :: candidate(:, :, :, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: rhs(:)
    integer :: i, j, k, previous_j, previous_k

    allocate(candidate, mold=old_state)
    allocate(rhs(size(old_state, 1)))
    do k = 1, size(old_state, 4)
      previous_k = modulo(k - 2, size(old_state, 4)) + 1
      do j = 1, size(old_state, 3)
        previous_j = modulo(j - 2, size(old_state, 3)) + 1
        do i = 1, size(old_state, 2)
          rhs = -(face_flux_x(:, i, j, k) - &
              face_flux_x(:, i - 1, j, k)) / dx &
            -(face_flux_y(:, i, j, k) - &
              face_flux_y(:, i, previous_j, k)) / dy &
            -(face_flux_z(:, i, j, k) - &
              face_flux_z(:, i, j, previous_k)) / dz
          if (second_stage) then
            candidate(:, i, j, k) = 0.5_dp * old_state(:, i, j, k) + &
              0.5_dp * (stage_base(:, i, j, k) + dt * rhs)
          else
            candidate(:, i, j, k) = &
              stage_base(:, i, j, k) + dt * rhs
          end if
        end do
      end do
    end do
    ok = all(ieee_is_finite(candidate))
  end subroutine update_sparse_periodic_stage_3d

  subroutine advance_sparse_fine_patch_ssprk2_3d( &
      distribution, species, patch, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, state, temperature, &
      dx, dy, dz, dt, alpha_start, alpha_end, riemann_solver, &
      reconstruction, limiter, face_flux_x, face_flux_y, face_flux_z, ok)
    type(mpi_amr_sparse_distribution_3d), intent(in) :: distribution
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: coarse_start(:, :, :, :)
    real(dp), intent(in) :: coarse_start_temperature(:, :, :)
    real(dp), intent(in) :: coarse_end(:, :, :, :)
    real(dp), intent(in) :: coarse_end_temperature(:, :, :)
    real(dp), intent(inout) :: state(:, :, :, :), temperature(:, :, :)
    real(dp), intent(in) :: dx, dy, dz, dt, alpha_start, alpha_end
    character(len=*), intent(in) :: riemann_solver, reconstruction, limiter
    real(dp), allocatable, intent(out) :: face_flux_x(:, :, :, :)
    real(dp), allocatable, intent(out) :: face_flux_y(:, :, :, :)
    real(dp), allocatable, intent(out) :: face_flux_z(:, :, :, :)
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
    logical :: local_ok
    integer :: nvar, local_count, ny, nz

    ok = .false.
    nvar = reactive_nvar(size(species))
    local_count = distribution%fine_count
    ny = patch%fine_ny()
    nz = patch%fine_nz()
    if (local_count == 0) then
      allocate(face_flux_x(nvar, 0:0, ny, nz))
      allocate(face_flux_y(nvar, 0, 0:ny, nz))
      allocate(face_flux_z(nvar, 0, ny, 0:nz))
      face_flux_x = 0.0_dp
      face_flux_y = 0.0_dp
      face_flux_z = 0.0_dp
      ok = size(state, 2) == 0 .and. size(temperature, 1) == 0
      return
    end if
    if (size(state, 1) /= nvar .or. size(state, 2) /= local_count .or. &
        size(state, 3) /= ny .or. size(state, 4) /= nz .or. &
        any(shape(temperature) /= [local_count, ny, nz])) return

    old_state = state
    allocate(old_temperature(local_count, ny, nz))
    call recover_reactive_temperatures_3d( &
      species, old_state, temperature, local_count, ny, nz, &
      old_temperature, local_ok)
    if (.not. local_ok) return
    call compute_sparse_fine_face_fluxes_3d( &
      distribution, species, patch, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, old_state, old_temperature, &
      alpha_start, riemann_solver, reconstruction, limiter, first_flux_x, &
      first_flux_y, first_flux_z, local_ok)
    if (.not. local_ok) return
    call update_sparse_fine_stage_3d( &
      old_state, old_state, first_flux_x, first_flux_y, first_flux_z, &
      dx, dy, dz, dt, .false., stage_state, local_ok)
    if (.not. local_ok) return
    allocate(stage_temperature(local_count, ny, nz))
    call recover_reactive_temperatures_3d( &
      species, stage_state, old_temperature, local_count, ny, nz, &
      stage_temperature, local_ok)
    if (.not. local_ok) return

    call compute_sparse_fine_face_fluxes_3d( &
      distribution, species, patch, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, stage_state, stage_temperature, &
      alpha_end, riemann_solver, reconstruction, limiter, second_flux_x, &
      second_flux_y, second_flux_z, local_ok)
    if (.not. local_ok) return
    call update_sparse_fine_stage_3d( &
      old_state, stage_state, second_flux_x, second_flux_y, second_flux_z, &
      dx, dy, dz, dt, .true., candidate_state, local_ok)
    if (.not. local_ok) return
    allocate(candidate_temperature(local_count, ny, nz))
    call recover_reactive_temperatures_3d( &
      species, candidate_state, stage_temperature, local_count, ny, nz, &
      candidate_temperature, local_ok)
    if (.not. local_ok) return

    allocate(face_flux_x(nvar, 0:local_count, ny, nz))
    allocate(face_flux_y(nvar, local_count, 0:ny, nz))
    allocate(face_flux_z(nvar, local_count, ny, 0:nz))
    face_flux_x = 0.5_dp * (first_flux_x + second_flux_x)
    face_flux_y = 0.5_dp * (first_flux_y + second_flux_y)
    face_flux_z = 0.5_dp * (first_flux_z + second_flux_z)
    state = candidate_state
    temperature = candidate_temperature
    ok = .true.
  end subroutine advance_sparse_fine_patch_ssprk2_3d

  subroutine compute_sparse_fine_face_fluxes_3d( &
      distribution, species, patch, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, state, temperature, alpha, &
      riemann_solver, reconstruction, limiter, face_flux_x, &
      face_flux_y, face_flux_z, ok)
    type(mpi_amr_sparse_distribution_3d), intent(in) :: distribution
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: coarse_start(:, :, :, :)
    real(dp), intent(in) :: coarse_start_temperature(:, :, :)
    real(dp), intent(in) :: coarse_end(:, :, :, :)
    real(dp), intent(in) :: coarse_end_temperature(:, :, :)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    real(dp), intent(in) :: alpha
    character(len=*), intent(in) :: riemann_solver, reconstruction, limiter
    real(dp), allocatable, intent(out) :: face_flux_x(:, :, :, :)
    real(dp), allocatable, intent(out) :: face_flux_y(:, :, :, :)
    real(dp), allocatable, intent(out) :: face_flux_z(:, :, :, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: extended_state(:, :, :, :)
    real(dp), allocatable :: extended_temperature(:, :, :)
    real(dp), allocatable :: extended_flux_x(:, :, :, :)
    real(dp), allocatable :: extended_flux_y(:, :, :, :)
    real(dp), allocatable :: extended_flux_z(:, :, :, :)
    logical :: face_ok
    integer :: i, j, k, nvar, local_count, ny, nz

    ok = .false.
    nvar = size(state, 1)
    local_count = size(state, 2)
    ny = size(state, 3)
    nz = size(state, 4)
    allocate(face_flux_x(nvar, 0:local_count, ny, nz))
    allocate(face_flux_y(nvar, local_count, 0:ny, nz))
    allocate(face_flux_z(nvar, local_count, ny, 0:nz))
    face_flux_x = 0.0_dp
    face_flux_y = 0.0_dp
    face_flux_z = 0.0_dp
    call build_sparse_fine_extended_state_3d( &
      distribution, species, patch, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, state, temperature, alpha, &
      reconstruction, limiter, extended_state, extended_temperature, face_ok)
    if (.not. face_ok) return

    select case (trim(reconstruction))
    case ("characteristic_plm")
      allocate(extended_flux_x, mold=extended_state)
      allocate(extended_flux_y, mold=extended_state)
      allocate(extended_flux_z, mold=extended_state)
      call compute_reactive_plm_slab_face_fluxes_3d( &
        species, extended_state, extended_temperature, &
        local_count + 4, ny + 4, nz + 4, 2, local_count + 2, &
        limiter, riemann_solver, extended_flux_x, extended_flux_y, &
        extended_flux_z, face_ok)
      if (.not. face_ok) return
      do k = 1, nz
        do j = 1, ny
          do i = 0, local_count
            face_flux_x(:, i, j, k) = &
              extended_flux_x(:, i + 2, j + 2, k + 2)
          end do
        end do
      end do
      do k = 1, nz
        do j = 0, ny
          do i = 1, local_count
            face_flux_y(:, i, j, k) = &
              extended_flux_y(:, i + 2, j + 2, k + 2)
          end do
        end do
      end do
      do k = 0, nz
        do j = 1, ny
          do i = 1, local_count
            face_flux_z(:, i, j, k) = &
              extended_flux_z(:, i + 2, j + 2, k + 2)
          end do
        end do
      end do
    case ("pcm")
      x_faces: do k = 1, nz
        do j = 1, ny
          do i = 0, local_count
            call reactive_riemann_flux_x( &
              species, extended_state(:, i + 2, j + 2, k + 2), &
              extended_state(:, i + 3, j + 2, k + 2), &
              extended_temperature(i + 2, j + 2, k + 2), &
              extended_temperature(i + 3, j + 2, k + 2), &
              riemann_solver, face_flux_x(:, i, j, k), face_ok)
            if (.not. face_ok) exit x_faces
          end do
        end do
      end do x_faces
      if (.not. face_ok) return
      y_faces: do k = 1, nz
        do j = 0, ny
          do i = 1, local_count
            call reactive_riemann_flux_y( &
              species, extended_state(:, i + 2, j + 2, k + 2), &
              extended_state(:, i + 2, j + 3, k + 2), &
              extended_temperature(i + 2, j + 2, k + 2), &
              extended_temperature(i + 2, j + 3, k + 2), &
              riemann_solver, face_flux_y(:, i, j, k), face_ok)
            if (.not. face_ok) exit y_faces
          end do
        end do
      end do y_faces
      if (.not. face_ok) return
      z_faces: do k = 0, nz
        do j = 1, ny
          do i = 1, local_count
            call reactive_riemann_flux_z( &
              species, extended_state(:, i + 2, j + 2, k + 2), &
              extended_state(:, i + 2, j + 2, k + 3), &
              extended_temperature(i + 2, j + 2, k + 2), &
              extended_temperature(i + 2, j + 2, k + 3), &
              riemann_solver, face_flux_z(:, i, j, k), face_ok)
            if (.not. face_ok) exit z_faces
          end do
        end do
      end do z_faces
      if (.not. face_ok) return
    case default
      return
    end select
    ok = all(ieee_is_finite(face_flux_x)) .and. &
      all(ieee_is_finite(face_flux_y)) .and. &
      all(ieee_is_finite(face_flux_z))
  end subroutine compute_sparse_fine_face_fluxes_3d

  subroutine build_sparse_fine_extended_state_3d( &
      distribution, species, patch, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, state, temperature, alpha, &
      reconstruction, limiter, extended_state, extended_temperature, ok)
    type(mpi_amr_sparse_distribution_3d), intent(in) :: distribution
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: coarse_start(:, :, :, :)
    real(dp), intent(in) :: coarse_start_temperature(:, :, :)
    real(dp), intent(in) :: coarse_end(:, :, :, :)
    real(dp), intent(in) :: coarse_end_temperature(:, :, :)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    real(dp), intent(in) :: alpha
    character(len=*), intent(in) :: reconstruction, limiter
    real(dp), allocatable, intent(out) :: extended_state(:, :, :, :)
    real(dp), allocatable, intent(out) :: extended_temperature(:, :, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: fine_x_halo_state(:, :, :, :)
    real(dp), allocatable :: fine_x_halo_temperature(:, :, :)
    real(dp), allocatable :: primitive(:)
    real(dp) :: ghost_temperature
    logical :: local_ok
    integer :: i, j, k, fine_i, fine_j, fine_k
    integer :: coarse_i, coarse_j, coarse_k
    integer :: local_count, ny, nz, ratio, nvar

    ok = .false.
    local_count = distribution%fine_count
    ny = patch%fine_ny()
    nz = patch%fine_nz()
    ratio = patch%refinement_ratio
    nvar = reactive_nvar(size(species))
    if (local_count < 1 .or. size(state, 1) /= nvar .or. &
        size(state, 2) /= local_count .or. size(state, 3) /= ny .or. &
        size(state, 4) /= nz) return
    if (.not. all(shape(temperature) == [local_count, ny, nz])) return
    if (.not. ieee_is_finite(alpha)) return
    if (alpha < 0.0_dp .or. alpha > 1.0_dp) return
    if (.not. all(ieee_is_finite(state))) return
    if (.not. all(ieee_is_finite(temperature))) return
    allocate(fine_x_halo_state(nvar, local_count + 4, ny, nz))
    allocate(fine_x_halo_temperature(local_count + 4, ny, nz))
    fine_x_halo_state = 0.0_dp
    fine_x_halo_temperature = 0.0_dp
    fine_x_halo_state(:, 3:local_count + 2, :, :) = state
    fine_x_halo_temperature(3:local_count + 2, :, :) = temperature
    call exchange_sparse_x_halos_3d( &
      distribution%comm, distribution%previous_fine_rank, &
      distribution%next_fine_rank, fine_x_halo_state, &
      fine_x_halo_temperature, local_count, local_ok, .false.)
    if (.not. local_ok) return

    allocate(extended_state(nvar, local_count + 4, ny + 4, nz + 4))
    allocate(extended_temperature(local_count + 4, ny + 4, nz + 4))
    allocate(primitive(reactive_nprim(size(species))))
    extended_state = 0.0_dp
    extended_temperature = 0.0_dp
    do k = 1, nz + 4
      fine_k = k - 2
      coarse_k = patch%coarse_k_lower + &
        floor(real(fine_k - 1, dp) / real(ratio, dp))
      do j = 1, ny + 4
        fine_j = j - 2
        coarse_j = patch%coarse_j_lower + &
          floor(real(fine_j - 1, dp) / real(ratio, dp))
        do i = 1, local_count + 4
          fine_i = distribution%fine_first + i - 3
          if (fine_i >= 1 .and. fine_i <= patch%fine_nx() .and. &
              fine_j >= 1 .and. fine_j <= ny .and. &
              fine_k >= 1 .and. fine_k <= nz) then
            extended_state(:, i, j, k) = fine_x_halo_state(:, i, fine_j, fine_k)
            extended_temperature(i, j, k) = &
              fine_x_halo_temperature(i, fine_j, fine_k)
            cycle
          end if
          coarse_i = patch%coarse_i_lower + &
            floor(real(fine_i - 1, dp) / real(ratio, dp))
          select case (trim(reconstruction))
          case ("characteristic_plm")
            call interpolate_sparse_coarse_fine_ghost_plm_3d( &
              distribution, species, patch, coarse_start, &
              coarse_start_temperature, coarse_end, coarse_end_temperature, &
              coarse_i, coarse_j, coarse_k, fine_i, fine_j, fine_k, &
              alpha, limiter, primitive, extended_state(:, i, j, k), &
              ghost_temperature, local_ok)
          case ("pcm")
            call interpolate_sparse_coarse_cell_3d( &
              distribution, species, patch, coarse_start, &
              coarse_start_temperature, coarse_end, coarse_end_temperature, &
              coarse_i, coarse_j, coarse_k, alpha, primitive, &
              extended_state(:, i, j, k), ghost_temperature, local_ok)
          case default
            local_ok = .false.
          end select
          if (.not. local_ok) return
          extended_temperature(i, j, k) = ghost_temperature
        end do
      end do
    end do
    ok = all(ieee_is_finite(extended_state))
    if (ok) ok = all(ieee_is_finite(extended_temperature))
    if (ok) ok = minval(extended_temperature) > 0.0_dp
  end subroutine build_sparse_fine_extended_state_3d

  subroutine update_sparse_fine_stage_3d( &
      old_state, stage_base, face_flux_x, face_flux_y, face_flux_z, &
      dx, dy, dz, dt, second_stage, candidate, ok)
    real(dp), intent(in) :: old_state(:, :, :, :)
    real(dp), intent(in) :: stage_base(:, :, :, :)
    real(dp), intent(in) :: face_flux_x(:, 0:, :, :)
    real(dp), intent(in) :: face_flux_y(:, :, 0:, :)
    real(dp), intent(in) :: face_flux_z(:, :, :, 0:)
    real(dp), intent(in) :: dx, dy, dz, dt
    logical, intent(in) :: second_stage
    real(dp), allocatable, intent(out) :: candidate(:, :, :, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: rhs(:)
    integer :: i, j, k

    allocate(candidate, mold=old_state)
    allocate(rhs(size(old_state, 1)))
    do k = 1, size(old_state, 4)
      do j = 1, size(old_state, 3)
        do i = 1, size(old_state, 2)
          rhs = -(face_flux_x(:, i, j, k) - &
              face_flux_x(:, i - 1, j, k)) / dx &
            -(face_flux_y(:, i, j, k) - &
              face_flux_y(:, i, j - 1, k)) / dy &
            -(face_flux_z(:, i, j, k) - &
              face_flux_z(:, i, j, k - 1)) / dz
          if (second_stage) then
            candidate(:, i, j, k) = 0.5_dp * old_state(:, i, j, k) + &
              0.5_dp * (stage_base(:, i, j, k) + dt * rhs)
          else
            candidate(:, i, j, k) = &
              stage_base(:, i, j, k) + dt * rhs
          end if
        end do
      end do
    end do
    ok = all(ieee_is_finite(candidate))
  end subroutine update_sparse_fine_stage_3d

  subroutine interpolate_sparse_coarse_cell_3d( &
      distribution, species, patch, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, coarse_i, coarse_j, coarse_k, &
      alpha, primitive, state, temperature, ok)
    type(mpi_amr_sparse_distribution_3d), intent(in) :: distribution
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: coarse_start(:, :, :, :)
    real(dp), intent(in) :: coarse_start_temperature(:, :, :)
    real(dp), intent(in) :: coarse_end(:, :, :, :)
    real(dp), intent(in) :: coarse_end_temperature(:, :, :)
    integer, intent(in) :: coarse_i, coarse_j, coarse_k
    real(dp), intent(in) :: alpha
    real(dp), intent(out) :: primitive(:), state(:), temperature
    logical, intent(out) :: ok

    real(dp) :: temperature_guess, sound_speed
    integer :: local_i, wrapped_j, wrapped_k

    local_i = sparse_halo_local_x_index_3d( &
      distribution, patch%coarse_nx, coarse_i)
    wrapped_j = 1 + modulo(coarse_j - 1, patch%coarse_ny)
    wrapped_k = 1 + modulo(coarse_k - 1, patch%coarse_nz)
    ok = local_i >= 1 .and. local_i <= size(coarse_start, 2)
    if (.not. ok) return
    state = (1.0_dp - alpha) * &
        coarse_start(:, local_i, wrapped_j, wrapped_k) + &
      alpha * coarse_end(:, local_i, wrapped_j, wrapped_k)
    temperature_guess = (1.0_dp - alpha) * &
        coarse_start_temperature(local_i, wrapped_j, wrapped_k) + &
      alpha * coarse_end_temperature(local_i, wrapped_j, wrapped_k)
    call reactive_conserved_to_primitive( &
      species, state, temperature_guess, primitive, temperature, &
      sound_speed, ok)
  end subroutine interpolate_sparse_coarse_cell_3d

  subroutine time_interpolated_sparse_coarse_primitive_3d( &
      distribution, species, patch, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, coarse_i, coarse_j, coarse_k, &
      alpha, primitive, ok)
    type(mpi_amr_sparse_distribution_3d), intent(in) :: distribution
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: coarse_start(:, :, :, :)
    real(dp), intent(in) :: coarse_start_temperature(:, :, :)
    real(dp), intent(in) :: coarse_end(:, :, :, :)
    real(dp), intent(in) :: coarse_end_temperature(:, :, :)
    integer, intent(in) :: coarse_i, coarse_j, coarse_k
    real(dp), intent(in) :: alpha
    real(dp), intent(out) :: primitive(:)
    logical, intent(out) :: ok

    real(dp) :: state(size(coarse_start, 1)), temperature

    call interpolate_sparse_coarse_cell_3d( &
      distribution, species, patch, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, coarse_i, coarse_j, coarse_k, &
      alpha, primitive, state, temperature, ok)
  end subroutine time_interpolated_sparse_coarse_primitive_3d

  subroutine interpolate_sparse_coarse_fine_ghost_plm_3d( &
      distribution, species, patch, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, coarse_i, coarse_j, coarse_k, &
      fine_i, fine_j, fine_k, alpha, limiter, primitive, state, &
      temperature, ok)
    type(mpi_amr_sparse_distribution_3d), intent(in) :: distribution
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: coarse_start(:, :, :, :)
    real(dp), intent(in) :: coarse_start_temperature(:, :, :)
    real(dp), intent(in) :: coarse_end(:, :, :, :)
    real(dp), intent(in) :: coarse_end_temperature(:, :, :)
    integer, intent(in) :: coarse_i, coarse_j, coarse_k
    integer, intent(in) :: fine_i, fine_j, fine_k
    real(dp), intent(in) :: alpha
    character(len=*), intent(in) :: limiter
    real(dp), intent(out) :: primitive(:), state(:), temperature
    logical, intent(out) :: ok

    real(dp) :: center(size(primitive)), minus(size(primitive))
    real(dp) :: plus(size(primitive)), slope_x(size(primitive))
    real(dp) :: slope_y(size(primitive)), slope_z(size(primitive))
    real(dp) :: delta(size(primitive))
    real(dp) :: offset_x, offset_y, offset_z, theta, candidate_value
    real(dp) :: sound_speed, mass_fraction_sum
    logical :: local_ok
    integer :: component, species_index, im, ip, jm, jp, km, kp, ratio

    ok = .false.
    ratio = patch%refinement_ratio
    im = 1 + modulo(coarse_i - 2, patch%coarse_nx)
    ip = 1 + modulo(coarse_i, patch%coarse_nx)
    jm = 1 + modulo(coarse_j - 2, patch%coarse_ny)
    jp = 1 + modulo(coarse_j, patch%coarse_ny)
    km = 1 + modulo(coarse_k - 2, patch%coarse_nz)
    kp = 1 + modulo(coarse_k, patch%coarse_nz)
    call time_interpolated_sparse_coarse_primitive_3d( &
      distribution, species, patch, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, coarse_i, coarse_j, coarse_k, &
      alpha, center, local_ok)
    if (.not. local_ok) return
    call time_interpolated_sparse_coarse_primitive_3d( &
      distribution, species, patch, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, im, coarse_j, coarse_k, &
      alpha, minus, local_ok)
    if (.not. local_ok) return
    call time_interpolated_sparse_coarse_primitive_3d( &
      distribution, species, patch, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, ip, coarse_j, coarse_k, &
      alpha, plus, local_ok)
    if (.not. local_ok) return
    do component = 1, size(primitive)
      call limited_slope( &
        center(component) - minus(component), &
        plus(component) - center(component), limiter, &
        slope_x(component), local_ok)
      if (.not. local_ok) return
    end do
    call time_interpolated_sparse_coarse_primitive_3d( &
      distribution, species, patch, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, coarse_i, jm, coarse_k, &
      alpha, minus, local_ok)
    if (.not. local_ok) return
    call time_interpolated_sparse_coarse_primitive_3d( &
      distribution, species, patch, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, coarse_i, jp, coarse_k, &
      alpha, plus, local_ok)
    if (.not. local_ok) return
    do component = 1, size(primitive)
      call limited_slope( &
        center(component) - minus(component), &
        plus(component) - center(component), limiter, &
        slope_y(component), local_ok)
      if (.not. local_ok) return
    end do
    call time_interpolated_sparse_coarse_primitive_3d( &
      distribution, species, patch, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, coarse_i, coarse_j, km, &
      alpha, minus, local_ok)
    if (.not. local_ok) return
    call time_interpolated_sparse_coarse_primitive_3d( &
      distribution, species, patch, coarse_start, coarse_start_temperature, &
      coarse_end, coarse_end_temperature, coarse_i, coarse_j, kp, &
      alpha, plus, local_ok)
    if (.not. local_ok) return
    do component = 1, size(primitive)
      call limited_slope( &
        center(component) - minus(component), &
        plus(component) - center(component), limiter, &
        slope_z(component), local_ok)
      if (.not. local_ok) return
    end do

    offset_x = (real(modulo(fine_i - 1, ratio), dp) + 0.5_dp) / &
      real(ratio, dp) - 0.5_dp
    offset_y = (real(modulo(fine_j - 1, ratio), dp) + 0.5_dp) / &
      real(ratio, dp) - 0.5_dp
    offset_z = (real(modulo(fine_k - 1, ratio), dp) + 0.5_dp) / &
      real(ratio, dp) - 0.5_dp
    delta = offset_x * slope_x + offset_y * slope_y + offset_z * slope_z
    theta = sparse_positive_increment_scale_3d( &
      center(1), delta(1), density_floor)
    theta = min(theta, sparse_positive_increment_scale_3d( &
      center(5), delta(5), pressure_floor))
    do species_index = 1, size(species)
      component = reactive_mass_fraction_component(species_index)
      theta = min(theta, sparse_positive_increment_scale_3d( &
        center(component), delta(component), 0.0_dp))
    end do
    primitive = center + max(0.0_dp, min(1.0_dp, theta)) * delta
    mass_fraction_sum = 0.0_dp
    do species_index = 1, size(species)
      component = reactive_mass_fraction_component(species_index)
      primitive(component) = max(0.0_dp, primitive(component))
      mass_fraction_sum = mass_fraction_sum + primitive(component)
    end do
    if (mass_fraction_sum <= tiny(1.0_dp)) return
    do species_index = 1, size(species)
      component = reactive_mass_fraction_component(species_index)
      primitive(component) = primitive(component) / mass_fraction_sum
    end do
    call reactive_primitive_to_conserved( &
      species, primitive, state, temperature, sound_speed, local_ok)
    if (.not. local_ok) then
      primitive = center
      call reactive_primitive_to_conserved( &
        species, primitive, state, temperature, sound_speed, local_ok)
    end if
    ok = local_ok
    if (ok) then
      candidate_value = temperature + sound_speed
      ok = ieee_is_finite(candidate_value)
    end if
  end subroutine interpolate_sparse_coarse_fine_ghost_plm_3d

  pure real(dp) function sparse_positive_increment_scale_3d( &
      center, delta, lower_bound) result(theta)
    real(dp), intent(in) :: center, delta, lower_bound

    theta = 1.0_dp
    if (delta < 0.0_dp) then
      theta = min(1.0_dp, max(0.0_dp, &
        (center - lower_bound) / max(-delta, tiny(1.0_dp))))
    end if
  end function sparse_positive_increment_scale_3d

  pure integer function sparse_halo_local_x_index_3d( &
      distribution, global_nx, global_i) result(local_i)
    type(mpi_amr_sparse_distribution_3d), intent(in) :: distribution
    integer, intent(in) :: global_nx, global_i

    integer :: wrapped_i, delta

    wrapped_i = 1 + modulo(global_i - 1, global_nx)
    delta = wrapped_i - distribution%coarse_first
    do while (delta < -2)
      delta = delta + global_nx
    end do
    do while (delta > distribution%coarse_count + 1)
      delta = delta - global_nx
    end do
    local_i = delta + 3
  end function sparse_halo_local_x_index_3d

  subroutine accumulate_sparse_fine_interface_fluxes_3d( &
      distribution, patch, face_flux_x, face_flux_y, face_flux_z, &
      x_lower, x_upper, y_lower, y_upper, z_lower, z_upper, ok)
    type(mpi_amr_sparse_distribution_3d), intent(in) :: distribution
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: face_flux_x(:, 0:, :, :)
    real(dp), intent(in) :: face_flux_y(:, :, 0:, :)
    real(dp), intent(in) :: face_flux_z(:, :, :, 0:)
    real(dp), intent(inout) :: x_lower(:, :, :), x_upper(:, :, :)
    real(dp), intent(inout) :: y_lower(:, :, :), y_upper(:, :, :)
    real(dp), intent(inout) :: z_lower(:, :, :), z_upper(:, :, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: local_x_lower(:, :, :), local_x_upper(:, :, :)
    real(dp), allocatable :: local_y_lower(:, :, :), local_y_upper(:, :, :)
    real(dp), allocatable :: local_z_lower(:, :, :), local_z_upper(:, :, :)
    real(dp), allocatable :: global_values(:, :, :)
    real(dp) :: inverse_face_children
    integer :: coarse_i, covered_i, fine_i, fine_j, fine_k
    integer :: j, k, local_fine_i, ratio, ierr
    logical :: collective_ok

    ok = .false.
    ratio = patch%refinement_ratio
    inverse_face_children = 1.0_dp / real(ratio**2, dp)
    allocate(local_x_lower, mold=x_lower)
    allocate(local_x_upper, mold=x_upper)
    allocate(local_y_lower, mold=y_lower)
    allocate(local_y_upper, mold=y_upper)
    allocate(local_z_lower, mold=z_lower)
    allocate(local_z_upper, mold=z_upper)
    local_x_lower = 0.0_dp
    local_x_upper = 0.0_dp
    local_y_lower = 0.0_dp
    local_y_upper = 0.0_dp
    local_z_lower = 0.0_dp
    local_z_upper = 0.0_dp

    if (distribution%fine_count > 0) then
      if (distribution%fine_first == 1) then
        local_x_lower = x_lower
        do k = 1, size(x_lower, 3)
          do j = 1, size(x_lower, 2)
            do fine_k = (k - 1) * ratio + 1, k * ratio
              do fine_j = (j - 1) * ratio + 1, j * ratio
                local_x_lower(:, j, k) = local_x_lower(:, j, k) + &
                  inverse_face_children * face_flux_x(:, 0, fine_j, fine_k)
              end do
            end do
          end do
        end do
      end if
      if (distribution%fine_last == patch%fine_nx()) then
        local_x_upper = x_upper
        do k = 1, size(x_upper, 3)
          do j = 1, size(x_upper, 2)
            do fine_k = (k - 1) * ratio + 1, k * ratio
              do fine_j = (j - 1) * ratio + 1, j * ratio
                local_x_upper(:, j, k) = local_x_upper(:, j, k) + &
                  inverse_face_children * face_flux_x( &
                    :, distribution%fine_count, fine_j, fine_k)
              end do
            end do
          end do
        end do
      end if

      do coarse_i = max(distribution%coarse_first, patch%coarse_i_lower), &
          min(distribution%coarse_last, patch%coarse_i_upper)
        covered_i = coarse_i - patch%coarse_i_lower + 1
        local_y_lower(:, covered_i, :) = y_lower(:, covered_i, :)
        local_y_upper(:, covered_i, :) = y_upper(:, covered_i, :)
        local_z_lower(:, covered_i, :) = z_lower(:, covered_i, :)
        local_z_upper(:, covered_i, :) = z_upper(:, covered_i, :)
        do k = 1, size(y_lower, 3)
          do fine_k = (k - 1) * ratio + 1, k * ratio
            do fine_i = (covered_i - 1) * ratio + 1, covered_i * ratio
              local_fine_i = fine_i - distribution%fine_first + 1
              local_y_lower(:, covered_i, k) = &
                local_y_lower(:, covered_i, k) + inverse_face_children * &
                  face_flux_y(:, local_fine_i, 0, fine_k)
              local_y_upper(:, covered_i, k) = &
                local_y_upper(:, covered_i, k) + inverse_face_children * &
                  face_flux_y(:, local_fine_i, patch%fine_ny(), fine_k)
            end do
          end do
        end do
        do j = 1, size(z_lower, 3)
          do fine_j = (j - 1) * ratio + 1, j * ratio
            do fine_i = (covered_i - 1) * ratio + 1, covered_i * ratio
              local_fine_i = fine_i - distribution%fine_first + 1
              local_z_lower(:, covered_i, j) = &
                local_z_lower(:, covered_i, j) + inverse_face_children * &
                  face_flux_z(:, local_fine_i, fine_j, 0)
              local_z_upper(:, covered_i, j) = &
                local_z_upper(:, covered_i, j) + inverse_face_children * &
                  face_flux_z(:, local_fine_i, fine_j, patch%fine_nz())
            end do
          end do
        end do
      end do
    end if

    allocate(global_values, mold=x_lower)
    call MPI_Allreduce(local_x_lower, global_values, size(x_lower), &
      MPI_DOUBLE_PRECISION, MPI_SUM, distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, collective_ok)
    if (.not. collective_ok) return
    x_lower = global_values
    call MPI_Allreduce(local_x_upper, global_values, size(x_upper), &
      MPI_DOUBLE_PRECISION, MPI_SUM, distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, collective_ok)
    if (.not. collective_ok) return
    x_upper = global_values
    deallocate(global_values)
    allocate(global_values, mold=y_lower)
    call MPI_Allreduce(local_y_lower, global_values, size(y_lower), &
      MPI_DOUBLE_PRECISION, MPI_SUM, distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, collective_ok)
    if (.not. collective_ok) return
    y_lower = global_values
    call MPI_Allreduce(local_y_upper, global_values, size(y_upper), &
      MPI_DOUBLE_PRECISION, MPI_SUM, distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, collective_ok)
    if (.not. collective_ok) return
    y_upper = global_values
    deallocate(global_values)
    allocate(global_values, mold=z_lower)
    call MPI_Allreduce(local_z_lower, global_values, size(z_lower), &
      MPI_DOUBLE_PRECISION, MPI_SUM, distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, collective_ok)
    if (.not. collective_ok) return
    z_lower = global_values
    call MPI_Allreduce(local_z_upper, global_values, size(z_upper), &
      MPI_DOUBLE_PRECISION, MPI_SUM, distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, collective_ok)
    if (.not. collective_ok) return
    z_upper = global_values
    ok = all(ieee_is_finite(x_lower)) .and. &
      all(ieee_is_finite(x_upper)) .and. &
      all(ieee_is_finite(y_lower)) .and. &
      all(ieee_is_finite(y_upper)) .and. &
      all(ieee_is_finite(z_lower)) .and. all(ieee_is_finite(z_upper))
  end subroutine accumulate_sparse_fine_interface_fluxes_3d

  subroutine reflux_sparse_coarse_3d( &
      distribution, patch, coarse_state, coarse_flux_x, coarse_flux_y, &
      coarse_flux_z, fine_x_lower, fine_x_upper, fine_y_lower, fine_y_upper, &
      fine_z_lower, fine_z_upper, dx, dy, dz, dt, maximum_correction, ok)
    type(mpi_amr_sparse_distribution_3d), intent(in) :: distribution
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(inout) :: coarse_state(:, :, :, :)
    real(dp), intent(in) :: coarse_flux_x(:, 0:, :, :)
    real(dp), intent(in) :: coarse_flux_y(:, :, :, :)
    real(dp), intent(in) :: coarse_flux_z(:, :, :, :)
    real(dp), intent(in) :: fine_x_lower(:, :, :), fine_x_upper(:, :, :)
    real(dp), intent(in) :: fine_y_lower(:, :, :), fine_y_upper(:, :, :)
    real(dp), intent(in) :: fine_z_lower(:, :, :), fine_z_upper(:, :, :)
    real(dp), intent(in) :: dx, dy, dz, dt
    real(dp), intent(out) :: maximum_correction
    logical, intent(out) :: ok

    real(dp), allocatable :: correction(:)
    integer :: global_i, local_i, local_j, local_k, j, k, face_i

    maximum_correction = 0.0_dp
    ok = .false.
    if (.not. all(ieee_is_finite([dx, dy, dz, dt]))) return
    if (dx <= 0.0_dp .or. dy <= 0.0_dp .or. dz <= 0.0_dp .or. &
        dt <= 0.0_dp) return
    if (.not. all(ieee_is_finite(coarse_state))) return
    if (.not. all(ieee_is_finite(coarse_flux_x))) return
    if (.not. all(ieee_is_finite(coarse_flux_y))) return
    if (.not. all(ieee_is_finite(coarse_flux_z))) return
    if (.not. all(ieee_is_finite(fine_x_lower))) return
    if (.not. all(ieee_is_finite(fine_x_upper))) return
    if (.not. all(ieee_is_finite(fine_y_lower))) return
    if (.not. all(ieee_is_finite(fine_y_upper))) return
    if (.not. all(ieee_is_finite(fine_z_lower))) return
    if (.not. all(ieee_is_finite(fine_z_upper))) return
    allocate(correction(size(coarse_state, 1)))
    global_i = patch%coarse_i_lower - 1
    if (distribution%owns_coarse(global_i)) then
      local_i = distribution%coarse_local_index(global_i)
      face_i = global_i - distribution%coarse_first + 1
      do k = patch%coarse_k_lower, patch%coarse_k_upper
        local_k = k - patch%coarse_k_lower + 1
        do j = patch%coarse_j_lower, patch%coarse_j_upper
          local_j = j - patch%coarse_j_lower + 1
          correction = -dt / dx * (fine_x_lower(:, local_j, local_k) - &
            coarse_flux_x(:, face_i, j, k))
          if (.not. all(ieee_is_finite(correction))) return
          coarse_state(:, local_i, j, k) = &
            coarse_state(:, local_i, j, k) + correction
          maximum_correction = max( &
            maximum_correction, maxval(abs(correction)))
        end do
      end do
    end if
    global_i = patch%coarse_i_upper + 1
    if (distribution%owns_coarse(global_i)) then
      local_i = distribution%coarse_local_index(global_i)
      face_i = patch%coarse_i_upper - distribution%coarse_first + 1
      do k = patch%coarse_k_lower, patch%coarse_k_upper
        local_k = k - patch%coarse_k_lower + 1
        do j = patch%coarse_j_lower, patch%coarse_j_upper
          local_j = j - patch%coarse_j_lower + 1
          correction = dt / dx * (fine_x_upper(:, local_j, local_k) - &
            coarse_flux_x(:, face_i, j, k))
          if (.not. all(ieee_is_finite(correction))) return
          coarse_state(:, local_i, j, k) = &
            coarse_state(:, local_i, j, k) + correction
          maximum_correction = max( &
            maximum_correction, maxval(abs(correction)))
        end do
      end do
    end if

    do global_i = max(distribution%coarse_first, patch%coarse_i_lower), &
        min(distribution%coarse_last, patch%coarse_i_upper)
      local_i = distribution%coarse_local_index(global_i)
      do k = patch%coarse_k_lower, patch%coarse_k_upper
        local_k = k - patch%coarse_k_lower + 1
        correction = -dt / dy * ( &
          fine_y_lower(:, global_i - patch%coarse_i_lower + 1, local_k) - &
            coarse_flux_y(:, local_i, patch%coarse_j_lower - 1, k))
        if (.not. all(ieee_is_finite(correction))) return
        coarse_state(:, local_i, patch%coarse_j_lower - 1, k) = &
          coarse_state(:, local_i, patch%coarse_j_lower - 1, k) + correction
        maximum_correction = max(maximum_correction, maxval(abs(correction)))
        correction = dt / dy * ( &
          fine_y_upper(:, global_i - patch%coarse_i_lower + 1, local_k) - &
            coarse_flux_y(:, local_i, patch%coarse_j_upper, k))
        if (.not. all(ieee_is_finite(correction))) return
        coarse_state(:, local_i, patch%coarse_j_upper + 1, k) = &
          coarse_state(:, local_i, patch%coarse_j_upper + 1, k) + correction
        maximum_correction = max(maximum_correction, maxval(abs(correction)))
      end do
      do j = patch%coarse_j_lower, patch%coarse_j_upper
        local_j = j - patch%coarse_j_lower + 1
        correction = -dt / dz * ( &
          fine_z_lower(:, global_i - patch%coarse_i_lower + 1, local_j) - &
            coarse_flux_z(:, local_i, j, patch%coarse_k_lower - 1))
        if (.not. all(ieee_is_finite(correction))) return
        coarse_state(:, local_i, j, patch%coarse_k_lower - 1) = &
          coarse_state(:, local_i, j, patch%coarse_k_lower - 1) + correction
        maximum_correction = max(maximum_correction, maxval(abs(correction)))
        correction = dt / dz * ( &
          fine_z_upper(:, global_i - patch%coarse_i_lower + 1, local_j) - &
            coarse_flux_z(:, local_i, j, patch%coarse_k_upper))
        if (.not. all(ieee_is_finite(correction))) return
        coarse_state(:, local_i, j, patch%coarse_k_upper + 1) = &
          coarse_state(:, local_i, j, patch%coarse_k_upper + 1) + correction
        maximum_correction = max(maximum_correction, maxval(abs(correction)))
      end do
    end do
    ok = ieee_is_finite(maximum_correction)
    if (ok) ok = all(ieee_is_finite(coarse_state))
  end subroutine reflux_sparse_coarse_3d

  subroutine average_down_sparse_hierarchy_3d( &
      distribution, patch, hierarchy, ok)
    type(mpi_amr_sparse_distribution_3d), intent(in) :: distribution
    type(amr_patch_3d), intent(in) :: patch
    type(mpi_amr_sparse_hierarchy_3d), intent(inout) :: hierarchy
    logical, intent(out) :: ok

    real(dp), allocatable :: restricted(:)
    real(dp) :: inverse_children
    integer :: global_i, local_i, covered_i, child_i, child_j, child_k
    integer :: local_child_i, j_lower, k_lower, coarse_j, coarse_k, ratio

    ok = .false.
    ratio = patch%refinement_ratio
    inverse_children = 1.0_dp / real(ratio**3, dp)
    allocate(restricted(hierarchy%nvar))
    do coarse_k = 1, patch%coarse_k_upper - patch%coarse_k_lower + 1
      k_lower = (coarse_k - 1) * ratio + 1
      do coarse_j = 1, patch%coarse_j_upper - patch%coarse_j_lower + 1
        j_lower = (coarse_j - 1) * ratio + 1
        do global_i = max( &
            distribution%coarse_first, patch%coarse_i_lower), &
            min(distribution%coarse_last, patch%coarse_i_upper)
          local_i = distribution%coarse_local_index(global_i)
          covered_i = global_i - patch%coarse_i_lower + 1
          restricted = 0.0_dp
          do child_k = k_lower, k_lower + ratio - 1
            do child_j = j_lower, j_lower + ratio - 1
              do child_i = (covered_i - 1) * ratio + 1, &
                  covered_i * ratio
                local_child_i = child_i - distribution%fine_first + 1
                restricted = restricted + inverse_children * &
                  hierarchy%fine_state( &
                    :, local_child_i, child_j, child_k)
              end do
            end do
          end do
          hierarchy%coarse_state(:, local_i, &
            patch%coarse_j_lower + coarse_j - 1, &
            patch%coarse_k_lower + coarse_k - 1) = restricted
        end do
      end do
    end do
    ok = all(ieee_is_finite(hierarchy%coarse_state))
  end subroutine average_down_sparse_hierarchy_3d

  subroutine compute_mpi_amr_sparse_transport_timestep_3d( &
      distribution, species, transport, patch, hierarchy, dx, dy, dz, &
      transport_cfl, viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, dt, maximum_diffusivity, ok)
    type(mpi_amr_sparse_distribution_3d), intent(in) :: distribution
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(amr_patch_3d), intent(in) :: patch
    type(mpi_amr_sparse_hierarchy_3d), intent(in) :: hierarchy
    real(dp), intent(in) :: dx, dy, dz, transport_cfl
    logical, intent(in) :: viscosity_enabled, thermal_conduction_enabled
    logical, intent(in) :: species_diffusion_enabled
    real(dp), intent(out) :: dt, maximum_diffusivity
    logical, intent(out) :: ok

    real(dp) :: local_dt(2), global_dt(2)
    real(dp) :: local_diffusivity(2), global_diffusivity(2), ratio_real
    logical :: local_ok, coarse_ok, fine_ok, global_ok
    integer :: ierr

    dt = 0.0_dp
    maximum_diffusivity = 0.0_dp
    ok = .false.
    local_ok = all(ieee_is_finite([dx, dy, dz, transport_cfl]))
    if (local_ok) local_ok = dx > 0.0_dp .and. dy > 0.0_dp .and. &
      dz > 0.0_dp .and. transport_cfl > 0.0_dp .and. &
      transport_cfl <= 0.5_dp
    if (local_ok) local_ok = hierarchy%nvar == reactive_nvar(size(species))
    if (local_ok) local_ok = hierarchy%is_valid(distribution, patch)
    if (local_ok) local_ok = compatible_transport_database(species, transport)
    call sparse_collective_contract_matches_3d( &
      distribution, species, patch, &
      [dx, dy, dz, transport_cfl, &
        merge(1.0_dp, 0.0_dp, viscosity_enabled), &
        merge(1.0_dp, 0.0_dp, thermal_conduction_enabled), &
        merge(1.0_dp, 0.0_dp, species_diffusion_enabled)], &
      "transport_timestep", local_ok, global_ok)
    if (.not. global_ok) return
    call sparse_collective_transport_matches_3d( &
      distribution, transport, .true., global_ok)
    if (.not. global_ok) return

    call reactive_transport_timestep_3d( &
      species, transport, hierarchy%coarse_state, &
      hierarchy%coarse_temperature, distribution%coarse_count, &
      patch%coarse_ny, patch%coarse_nz, dx, dy, dz, transport_cfl, &
      viscosity_enabled, thermal_conduction_enabled, &
      species_diffusion_enabled, local_dt(1), local_diffusivity(1), coarse_ok)
    ratio_real = real(patch%refinement_ratio, dp)
    if (distribution%fine_count > 0) then
      call reactive_transport_timestep_3d( &
        species, transport, hierarchy%fine_state, hierarchy%fine_temperature, &
        distribution%fine_count, patch%fine_ny(), patch%fine_nz(), &
        dx / ratio_real, dy / ratio_real, dz / ratio_real, transport_cfl, &
        viscosity_enabled, thermal_conduction_enabled, &
        species_diffusion_enabled, local_dt(2), local_diffusivity(2), fine_ok)
    else
      local_dt(2) = huge(1.0_dp)
      local_diffusivity(2) = 0.0_dp
      fine_ok = .true.
    end if
    local_ok = coarse_ok .and. fine_ok
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    call MPI_Allreduce( &
      local_dt, global_dt, size(local_dt), MPI_DOUBLE_PRECISION, MPI_MIN, &
      distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, global_ok)
    if (.not. global_ok) return
    call MPI_Allreduce( &
      local_diffusivity, global_diffusivity, size(local_diffusivity), &
      MPI_DOUBLE_PRECISION, MPI_MAX, distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, global_ok)
    if (.not. global_ok) return

    dt = min(global_dt(1), ratio_real * ratio_real * global_dt(2))
    maximum_diffusivity = max(global_diffusivity(1), global_diffusivity(2))
    local_ok = all(ieee_is_finite([dt, maximum_diffusivity]))
    if (local_ok) local_ok = dt > 0.0_dp
    if (local_ok) local_ok = maximum_diffusivity >= 0.0_dp
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    ok = global_ok
    if (.not. ok) then
      dt = 0.0_dp
      maximum_diffusivity = 0.0_dp
    end if
  end subroutine compute_mpi_amr_sparse_transport_timestep_3d

  subroutine compute_mpi_amr_sparse_cfl_timestep_3d( &
      distribution, species, patch, hierarchy, dx, dy, dz, cfl, dt, ok)
    type(mpi_amr_sparse_distribution_3d), intent(in) :: distribution
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    type(mpi_amr_sparse_hierarchy_3d), intent(in) :: hierarchy
    real(dp), intent(in) :: dx, dy, dz, cfl
    real(dp), intent(out) :: dt
    logical, intent(out) :: ok

    real(dp) :: local_coarse_rate, local_fine_rate
    real(dp) :: coarse_rate, fine_rate, ratio_real
    logical :: local_ok, local_coarse_ok, local_fine_ok, global_ok
    integer :: ierr

    dt = 0.0_dp
    ok = .false.
    local_ok = all(ieee_is_finite([dx, dy, dz, cfl]))
    if (local_ok) local_ok = dx > 0.0_dp .and. dy > 0.0_dp .and. &
      dz > 0.0_dp .and. cfl > 0.0_dp .and. cfl <= 1.0_dp
    if (local_ok) local_ok = hierarchy%nvar == reactive_nvar(size(species))
    if (local_ok) local_ok = hierarchy%is_valid(distribution, patch)
    call sparse_collective_contract_matches_3d( &
      distribution, species, patch, [dx, dy, dz, cfl], "", &
      local_ok, global_ok)
    if (.not. global_ok) return
    call local_sparse_maximum_rate_3d( &
      species, hierarchy%coarse_state, hierarchy%coarse_temperature, &
      dx, dy, dz, local_coarse_rate, local_coarse_ok)
    if (distribution%fine_count > 0) then
      ratio_real = real(patch%refinement_ratio, dp)
      call local_sparse_maximum_rate_3d( &
        species, hierarchy%fine_state, hierarchy%fine_temperature, &
        dx / ratio_real, dy / ratio_real, dz / ratio_real, &
        local_fine_rate, local_fine_ok)
    else
      local_fine_rate = 0.0_dp
      local_fine_ok = .true.
    end if
    local_ok = local_coarse_ok .and. local_fine_ok
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    call MPI_Allreduce( &
      local_coarse_rate, coarse_rate, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
      distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, global_ok)
    if (.not. global_ok) return
    call MPI_Allreduce( &
      local_fine_rate, fine_rate, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
      distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, global_ok)
    if (.not. global_ok) return
    ratio_real = real(patch%refinement_ratio, dp)
    local_ok = all(ieee_is_finite([coarse_rate, fine_rate]))
    if (local_ok) local_ok = coarse_rate > 0.0_dp .and. fine_rate > 0.0_dp
    if (local_ok) then
      dt = min(cfl / coarse_rate, ratio_real * cfl / fine_rate)
      local_ok = ieee_is_finite(dt)
      if (local_ok) local_ok = dt > 0.0_dp
    end if
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    ok = global_ok
    if (.not. ok) dt = 0.0_dp
  end subroutine compute_mpi_amr_sparse_cfl_timestep_3d

  subroutine local_sparse_maximum_rate_3d( &
      species, state, temperature, dx, dy, dz, maximum_rate, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    real(dp), intent(in) :: dx, dy, dz
    real(dp), intent(out) :: maximum_rate
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:)
    real(dp) :: recovered_temperature, sound_speed, rate
    logical :: cell_ok
    integer :: i, j, k

    maximum_rate = 0.0_dp
    ok = size(state, 2) > 0 .and. &
      all(shape(temperature) == [ &
        size(state, 2), size(state, 3), size(state, 4)])
    if (.not. ok) return
    ok = all(ieee_is_finite([dx, dy, dz]))
    if (.not. ok) return
    if (dx <= 0.0_dp .or. dy <= 0.0_dp .or. dz <= 0.0_dp) return
    ok = all(ieee_is_finite(state))
    if (.not. ok) return
    ok = all(ieee_is_finite(temperature))
    if (.not. ok) return
    allocate(primitive(reactive_nprim(size(species))))
    outer: do k = 1, size(state, 4)
      do j = 1, size(state, 3)
        do i = 1, size(state, 2)
          call reactive_conserved_to_primitive( &
            species, state(:, i, j, k), temperature(i, j, k), primitive, &
            recovered_temperature, sound_speed, cell_ok)
          if (.not. cell_ok) then
            ok = .false.
            exit outer
          end if
          rate = (abs(primitive(2)) + sound_speed) / dx + &
            (abs(primitive(3)) + sound_speed) / dy + &
            (abs(primitive(4)) + sound_speed) / dz
          if (.not. ieee_is_finite(rate)) then
            ok = .false.
            exit outer
          end if
          maximum_rate = max(maximum_rate, rate)
        end do
      end do
    end do outer
    if (.not. ok) return
    ok = ieee_is_finite(maximum_rate)
    if (ok) ok = maximum_rate > 0.0_dp
  end subroutine local_sparse_maximum_rate_3d

  subroutine build_mpi_amr_sparse_periodic_halos_3d( &
      distribution, state, temperature, halo_state, halo_temperature, ok)
    type(mpi_amr_sparse_distribution_3d), intent(in) :: distribution
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    real(dp), allocatable, intent(out) :: halo_state(:, :, :, :)
    real(dp), allocatable, intent(out) :: halo_temperature(:, :, :)
    logical, intent(out) :: ok

    logical :: local_ok, global_ok
    integer :: local_count, left_rank, right_rank, ierr
    integer :: local_shape(3), minimum_shape(3), maximum_shape(3)

    ok = .false.
    local_count = size(state, 2)
    local_shape = [size(state, 1), size(state, 3), size(state, 4)]
    call MPI_Allreduce( &
      local_shape, minimum_shape, size(local_shape), MPI_INTEGER, MPI_MIN, &
      distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, global_ok)
    if (.not. global_ok) return
    call MPI_Allreduce( &
      local_shape, maximum_shape, size(local_shape), MPI_INTEGER, MPI_MAX, &
      distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, global_ok)
    if (.not. global_ok) return
    local_ok = local_count >= 1 .and. &
      all(local_shape == minimum_shape) .and. &
      all(local_shape == maximum_shape) .and. &
      all(shape(temperature) == [ &
        local_count, size(state, 3), size(state, 4)])
    if (local_ok) local_ok = all(ieee_is_finite(state))
    if (local_ok) local_ok = all(ieee_is_finite(temperature))
    if (local_ok) local_ok = minval(temperature) > 0.0_dp
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    allocate(halo_state( &
      size(state, 1), local_count + 2 * mpi_amr_sparse_ghost_width_3d, &
      size(state, 3), size(state, 4)))
    allocate(halo_temperature( &
      local_count + 2 * mpi_amr_sparse_ghost_width_3d, &
      size(state, 3), size(state, 4)))
    halo_state = 0.0_dp
    halo_temperature = 0.0_dp
    halo_state(:, 3:local_count + 2, :, :) = state
    halo_temperature(3:local_count + 2, :, :) = temperature
    left_rank = modulo(distribution%rank - 1, distribution%nranks)
    right_rank = modulo(distribution%rank + 1, distribution%nranks)
    call exchange_sparse_x_halos_3d( &
      distribution%comm, left_rank, right_rank, halo_state, &
      halo_temperature, local_count, local_ok)
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    ok = global_ok
  end subroutine build_mpi_amr_sparse_periodic_halos_3d

  subroutine exchange_sparse_x_halos_3d( &
      comm, left_rank, right_rank, halo_state, halo_temperature, &
      local_count, ok, require_complete)
    type(MPI_Comm), intent(in) :: comm
    integer, intent(in) :: left_rank, right_rank, local_count
    real(dp), intent(inout) :: halo_state(:, :, :, :)
    real(dp), intent(inout) :: halo_temperature(:, :, :)
    logical, intent(out) :: ok
    logical, intent(in), optional :: require_complete

    real(dp), allocatable :: send_values(:), receive_values(:)
    integer :: depth, payload_size, send_index, receive_index, ierr
    type(MPI_Status) :: status
    logical :: complete

    ok = .false.
    complete = .true.
    if (present(require_complete)) complete = require_complete
    payload_size = (size(halo_state, 1) + 1) * &
      size(halo_state, 3) * size(halo_state, 4)
    allocate(send_values(payload_size), receive_values(payload_size))
    do depth = 1, mpi_amr_sparse_ghost_width_3d
      send_index = 3 + depth - 1
      receive_index = local_count + 2 + depth
      call pack_sparse_x_plane_3d( &
        halo_state, halo_temperature, send_index, send_values)
      receive_values = 0.0_dp
      call MPI_Sendrecv( &
        send_values, payload_size, MPI_DOUBLE_PRECISION, left_rank, &
        sparse_halo_left_tag_3d + depth, receive_values, payload_size, &
        MPI_DOUBLE_PRECISION, right_rank, sparse_halo_left_tag_3d + depth, &
        comm, status, ierr)
      if (ierr /= MPI_SUCCESS) return
      call unpack_sparse_x_plane_3d( &
        receive_values, halo_state, halo_temperature, receive_index)

      send_index = local_count + 2 - depth + 1
      receive_index = 3 - depth
      call pack_sparse_x_plane_3d( &
        halo_state, halo_temperature, send_index, send_values)
      receive_values = 0.0_dp
      call MPI_Sendrecv( &
        send_values, payload_size, MPI_DOUBLE_PRECISION, right_rank, &
        sparse_halo_right_tag_3d + depth, receive_values, payload_size, &
        MPI_DOUBLE_PRECISION, left_rank, sparse_halo_right_tag_3d + depth, &
        comm, status, ierr)
      if (ierr /= MPI_SUCCESS) return
      call unpack_sparse_x_plane_3d( &
        receive_values, halo_state, halo_temperature, receive_index)
    end do
    ok = all(ieee_is_finite(halo_state))
    if (ok) ok = all(ieee_is_finite(halo_temperature))
    if (complete) then
      if (ok) ok = minval(halo_temperature) > 0.0_dp
    end if
  end subroutine exchange_sparse_x_halos_3d

  subroutine pack_sparse_x_plane_3d( &
      state, temperature, plane, values)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: plane
    real(dp), intent(out) :: values(:)

    integer :: state_count

    state_count = size(state, 1) * size(state, 3) * size(state, 4)
    values(:state_count) = reshape( &
      state(:, plane, :, :), [state_count])
    values(state_count + 1:) = reshape( &
      temperature(plane, :, :), [size(temperature, 2) * &
        size(temperature, 3)])
  end subroutine pack_sparse_x_plane_3d

  subroutine unpack_sparse_x_plane_3d( &
      values, state, temperature, plane)
    real(dp), intent(in) :: values(:)
    real(dp), intent(inout) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: plane

    integer :: state_count

    state_count = size(state, 1) * size(state, 3) * size(state, 4)
    state(:, plane, :, :) = reshape( &
      values(:state_count), &
      [size(state, 1), size(state, 3), size(state, 4)])
    temperature(plane, :, :) = reshape( &
      values(state_count + 1:), &
      [size(temperature, 2), size(temperature, 3)])
  end subroutine unpack_sparse_x_plane_3d

  subroutine initialize_mpi_amr_sparse_distribution_3d( &
      patch, comm, distribution, ok)
    type(amr_patch_3d), intent(in) :: patch
    type(MPI_Comm), intent(in) :: comm
    type(mpi_amr_sparse_distribution_3d), intent(out) :: distribution
    logical, intent(out) :: ok

    integer :: ierr, rank_index, candidate_rank, parent_first, parent_last

    ok = .false.
    distribution%comm = comm
    call MPI_Comm_rank(comm, distribution%rank, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Comm_size(comm, distribution%nranks, ierr)
    if (ierr /= MPI_SUCCESS) return
    if (.not. patch%is_strictly_interior() .or. &
        distribution%nranks < 1 .or. &
        distribution%nranks > patch%coarse_nx) return

    allocate(distribution%coarse_counts(distribution%nranks))
    allocate(distribution%coarse_displacements(distribution%nranks))
    allocate(distribution%fine_counts(distribution%nranks))
    allocate(distribution%fine_displacements(distribution%nranks))
    call build_sparse_slab_counts_3d( &
      patch%coarse_nx, distribution%nranks, &
      distribution%coarse_counts, distribution%coarse_displacements)
    distribution%fine_counts = 0
    distribution%fine_displacements = 0
    do rank_index = 1, distribution%nranks
      parent_first = max( &
        distribution%coarse_displacements(rank_index) + 1, &
        patch%coarse_i_lower)
      parent_last = min( &
        distribution%coarse_displacements(rank_index) + &
          distribution%coarse_counts(rank_index), &
        patch%coarse_i_upper)
      if (parent_last >= parent_first) then
        distribution%fine_counts(rank_index) = &
          (parent_last - parent_first + 1) * patch%refinement_ratio
      end if
      if (rank_index > 1) then
        distribution%fine_displacements(rank_index) = &
          distribution%fine_displacements(rank_index - 1) + &
          distribution%fine_counts(rank_index - 1)
      end if
    end do

    rank_index = distribution%rank + 1
    distribution%coarse_count = distribution%coarse_counts(rank_index)
    distribution%coarse_first = &
      distribution%coarse_displacements(rank_index) + 1
    distribution%coarse_last = distribution%coarse_first + &
      distribution%coarse_count - 1
    distribution%fine_count = distribution%fine_counts(rank_index)
    distribution%previous_fine_rank = MPI_PROC_NULL
    distribution%next_fine_rank = MPI_PROC_NULL
    if (distribution%fine_count > 0) then
      distribution%fine_first = &
        distribution%fine_displacements(rank_index) + 1
      distribution%fine_last = distribution%fine_first + &
        distribution%fine_count - 1
      do candidate_rank = distribution%rank - 1, 0, -1
        if (distribution%fine_counts(candidate_rank + 1) > 0) then
          distribution%previous_fine_rank = candidate_rank
          exit
        end if
      end do
      do candidate_rank = distribution%rank + 1, distribution%nranks - 1
        if (distribution%fine_counts(candidate_rank + 1) > 0) then
          distribution%next_fine_rank = candidate_rank
          exit
        end if
      end do
    end if
    ok = distribution%is_valid(patch)
  end subroutine initialize_mpi_amr_sparse_distribution_3d

  pure subroutine build_sparse_slab_counts_3d( &
      global_count, nranks, counts, displacements)
    integer, intent(in) :: global_count, nranks
    integer, intent(out) :: counts(:), displacements(:)

    integer :: rank_index, base_count, remainder

    base_count = global_count / nranks
    remainder = modulo(global_count, nranks)
    counts = base_count
    do rank_index = 1, remainder
      counts(rank_index) = counts(rank_index) + 1
    end do
    displacements(1) = 0
    do rank_index = 2, nranks
      displacements(rank_index) = displacements(rank_index - 1) + &
        counts(rank_index - 1)
    end do
  end subroutine build_sparse_slab_counts_3d

  pure logical function sparse_distribution_is_valid_3d( &
      self, patch) result(valid)
    class(mpi_amr_sparse_distribution_3d), intent(in) :: self
    type(amr_patch_3d), intent(in) :: patch

    integer :: rank_index, candidate_rank, expected_fine_count
    integer :: expected_previous_fine_rank, expected_next_fine_rank

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
      all(self%coarse_counts >= 1) .and. all(self%fine_counts >= 0) .and. &
      sum(self%coarse_counts) == patch%coarse_nx .and. &
      sum(self%fine_counts) == patch%fine_nx() .and. &
      self%coarse_displacements(1) == 0 .and. &
      self%fine_displacements(1) == 0
    if (.not. valid) return
    do rank_index = 2, self%nranks
      if (self%coarse_displacements(rank_index) /= &
          self%coarse_displacements(rank_index - 1) + &
            self%coarse_counts(rank_index - 1)) then
        valid = .false.
        return
      end if
      if (self%fine_displacements(rank_index) /= &
          self%fine_displacements(rank_index - 1) + &
            self%fine_counts(rank_index - 1)) then
        valid = .false.
        return
      end if
    end do
    rank_index = self%rank + 1
    valid = self%coarse_count == self%coarse_counts(rank_index) .and. &
      self%coarse_first == self%coarse_displacements(rank_index) + 1 .and. &
      self%coarse_last == self%coarse_first + self%coarse_count - 1 .and. &
      self%fine_count == self%fine_counts(rank_index)
    if (.not. valid) return
    expected_fine_count = max(0, &
      min(self%coarse_last, patch%coarse_i_upper) - &
        max(self%coarse_first, patch%coarse_i_lower) + 1) * &
      patch%refinement_ratio
    valid = self%fine_count == expected_fine_count
    if (.not. valid) return
    if (self%fine_count > 0) then
      expected_previous_fine_rank = MPI_PROC_NULL
      expected_next_fine_rank = MPI_PROC_NULL
      do candidate_rank = self%rank - 1, 0, -1
        if (self%fine_counts(candidate_rank + 1) > 0) then
          expected_previous_fine_rank = candidate_rank
          exit
        end if
      end do
      do candidate_rank = self%rank + 1, self%nranks - 1
        if (self%fine_counts(candidate_rank + 1) > 0) then
          expected_next_fine_rank = candidate_rank
          exit
        end if
      end do
      valid = self%fine_first == self%fine_displacements(rank_index) + 1 .and. &
        self%fine_last == self%fine_first + self%fine_count - 1 .and. &
        self%fine_first >= 1 .and. self%fine_last <= patch%fine_nx() .and. &
        self%previous_fine_rank == expected_previous_fine_rank .and. &
        self%next_fine_rank == expected_next_fine_rank
    else
      valid = self%fine_first == 0 .and. self%fine_last == -1 .and. &
        self%previous_fine_rank == MPI_PROC_NULL .and. &
        self%next_fine_rank == MPI_PROC_NULL
    end if
  end function sparse_distribution_is_valid_3d

  pure logical function sparse_distribution_owns_coarse_3d( &
      self, global_i) result(owns)
    class(mpi_amr_sparse_distribution_3d), intent(in) :: self
    integer, intent(in) :: global_i

    owns = global_i >= self%coarse_first .and. &
      global_i <= self%coarse_last
  end function sparse_distribution_owns_coarse_3d

  pure logical function sparse_distribution_owns_fine_3d( &
      self, global_i) result(owns)
    class(mpi_amr_sparse_distribution_3d), intent(in) :: self
    integer, intent(in) :: global_i

    owns = self%fine_count > 0 .and. global_i >= self%fine_first .and. &
      global_i <= self%fine_last
  end function sparse_distribution_owns_fine_3d

  pure integer function sparse_coarse_local_index_3d( &
      self, global_i) result(local_i)
    class(mpi_amr_sparse_distribution_3d), intent(in) :: self
    integer, intent(in) :: global_i

    local_i = global_i - self%coarse_first + 1
  end function sparse_coarse_local_index_3d

  pure integer function sparse_fine_local_index_3d( &
      self, global_i) result(local_i)
    class(mpi_amr_sparse_distribution_3d), intent(in) :: self
    integer, intent(in) :: global_i

    local_i = global_i - self%fine_first + 1
  end function sparse_fine_local_index_3d

  pure logical function sparse_hierarchy_is_valid_3d( &
      self, distribution, patch) result(valid)
    class(mpi_amr_sparse_hierarchy_3d), intent(in) :: self
    type(mpi_amr_sparse_distribution_3d), intent(in) :: distribution
    type(amr_patch_3d), intent(in) :: patch

    valid = distribution%is_valid(patch) .and. self%nvar >= 6 .and. &
      allocated(self%coarse_state) .and. &
      allocated(self%coarse_temperature) .and. &
      allocated(self%fine_state) .and. allocated(self%fine_temperature)
    if (.not. valid) return
    valid = all(shape(self%coarse_state) == [ &
        self%nvar, distribution%coarse_count, &
        patch%coarse_ny, patch%coarse_nz]) .and. &
      all(shape(self%coarse_temperature) == [ &
        distribution%coarse_count, patch%coarse_ny, patch%coarse_nz]) .and. &
      all(shape(self%fine_state) == [ &
        self%nvar, distribution%fine_count, &
        patch%fine_ny(), patch%fine_nz()]) .and. &
      all(shape(self%fine_temperature) == [ &
        distribution%fine_count, patch%fine_ny(), patch%fine_nz()])
    if (.not. valid) return
    valid = all(ieee_is_finite(self%coarse_state))
    if (.not. valid) return
    valid = all(ieee_is_finite(self%coarse_temperature))
    if (.not. valid) return
    valid = all(ieee_is_finite(self%fine_state))
    if (.not. valid) return
    valid = all(ieee_is_finite(self%fine_temperature))
    if (.not. valid) return
    valid = minval(self%coarse_temperature) > 0.0_dp
    if (distribution%fine_count > 0) then
      if (valid) valid = minval(self%fine_temperature) > 0.0_dp
    end if
  end function sparse_hierarchy_is_valid_3d

  pure integer function sparse_hierarchy_local_cell_count_3d( &
      self) result(count)
    class(mpi_amr_sparse_hierarchy_3d), intent(in) :: self

    count = 0
    if (allocated(self%coarse_temperature)) &
      count = count + size(self%coarse_temperature)
    if (allocated(self%fine_temperature)) &
      count = count + size(self%fine_temperature)
  end function sparse_hierarchy_local_cell_count_3d

  pure integer function sparse_hierarchy_local_value_count_3d( &
      self) result(count)
    class(mpi_amr_sparse_hierarchy_3d), intent(in) :: self

    count = 0
    if (allocated(self%coarse_state)) count = count + size(self%coarse_state)
    if (allocated(self%coarse_temperature)) &
      count = count + size(self%coarse_temperature)
    if (allocated(self%fine_state)) count = count + size(self%fine_state)
    if (allocated(self%fine_temperature)) &
      count = count + size(self%fine_temperature)
  end function sparse_hierarchy_local_value_count_3d

  subroutine scatter_mpi_amr_sparse_hierarchy_3d( &
      distribution, species, patch, root, global_coarse_state, &
      global_coarse_temperature, global_fine_state, &
      global_fine_temperature, hierarchy, ok)
    type(mpi_amr_sparse_distribution_3d), intent(in) :: distribution
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    integer, intent(in) :: root
    real(dp), intent(in) :: global_coarse_state(:, :, :, :)
    real(dp), intent(in) :: global_coarse_temperature(:, :, :)
    real(dp), intent(in) :: global_fine_state(:, :, :, :)
    real(dp), intent(in) :: global_fine_temperature(:, :, :)
    type(mpi_amr_sparse_hierarchy_3d), intent(out) :: hierarchy
    logical, intent(out) :: ok

    type(mpi_amr_sparse_hierarchy_3d) :: candidate
    logical :: local_ok, global_ok
    integer :: nvar

    ok = .false.
    nvar = reactive_nvar(size(species))
    local_ok = distribution%is_valid(patch) .and. root >= 0 .and. &
      root < distribution%nranks .and. nvar >= 6
    call sparse_collective_root_matches_3d( &
      distribution, root, local_ok, global_ok)
    if (.not. global_ok) return
    call sparse_collective_contract_matches_3d( &
      distribution, species, patch, [real(nvar, dp)], "scatter", &
      local_ok, global_ok)
    if (.not. global_ok) return
    if (distribution%rank == root) then
      local_ok = local_ok .and. all(shape(global_coarse_state) == [ &
          nvar, patch%coarse_nx, patch%coarse_ny, patch%coarse_nz]) .and. &
        all(shape(global_coarse_temperature) == [ &
          patch%coarse_nx, patch%coarse_ny, patch%coarse_nz]) .and. &
        all(shape(global_fine_state) == [ &
          nvar, patch%fine_nx(), patch%fine_ny(), patch%fine_nz()]) .and. &
        all(shape(global_fine_temperature) == [ &
          patch%fine_nx(), patch%fine_ny(), patch%fine_nz()]) .and. &
        all(ieee_is_finite(global_coarse_state)) .and. &
        all(ieee_is_finite(global_coarse_temperature)) .and. &
        all(ieee_is_finite(global_fine_state)) .and. &
        all(ieee_is_finite(global_fine_temperature))
    end if
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return

    candidate%nvar = nvar
    allocate(candidate%coarse_state( &
      nvar, distribution%coarse_count, patch%coarse_ny, patch%coarse_nz))
    allocate(candidate%coarse_temperature( &
      distribution%coarse_count, patch%coarse_ny, patch%coarse_nz))
    allocate(candidate%fine_state( &
      nvar, distribution%fine_count, patch%fine_ny(), patch%fine_nz()))
    allocate(candidate%fine_temperature( &
      distribution%fine_count, patch%fine_ny(), patch%fine_nz()))
    call scatter_rank4_sparse_slabs_3d( &
      distribution%comm, distribution%rank, root, &
      distribution%coarse_counts, distribution%coarse_displacements, &
      global_coarse_state, candidate%coarse_state, local_ok)
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    call scatter_rank3_sparse_slabs_3d( &
      distribution%comm, distribution%rank, root, &
      distribution%coarse_counts, distribution%coarse_displacements, &
      global_coarse_temperature, candidate%coarse_temperature, local_ok)
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    call scatter_rank4_sparse_slabs_3d( &
      distribution%comm, distribution%rank, root, &
      distribution%fine_counts, distribution%fine_displacements, &
      global_fine_state, candidate%fine_state, local_ok)
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    call scatter_rank3_sparse_slabs_3d( &
      distribution%comm, distribution%rank, root, &
      distribution%fine_counts, distribution%fine_displacements, &
      global_fine_temperature, candidate%fine_temperature, local_ok)
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    local_ok = candidate%is_valid(distribution, patch)
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    hierarchy = candidate
    ok = .true.
  end subroutine scatter_mpi_amr_sparse_hierarchy_3d

  subroutine gather_mpi_amr_sparse_hierarchy_3d( &
      distribution, patch, root, hierarchy, global_coarse_state, &
      global_coarse_temperature, global_fine_state, &
      global_fine_temperature, ok)
    type(mpi_amr_sparse_distribution_3d), intent(in) :: distribution
    type(amr_patch_3d), intent(in) :: patch
    integer, intent(in) :: root
    type(mpi_amr_sparse_hierarchy_3d), intent(in) :: hierarchy
    real(dp), allocatable, intent(out) :: global_coarse_state(:, :, :, :)
    real(dp), allocatable, intent(out) :: global_coarse_temperature(:, :, :)
    real(dp), allocatable, intent(out) :: global_fine_state(:, :, :, :)
    real(dp), allocatable, intent(out) :: global_fine_temperature(:, :, :)
    logical, intent(out) :: ok

    logical :: local_ok, global_ok

    ok = .false.
    local_ok = root >= 0 .and. root < distribution%nranks .and. &
      hierarchy%is_valid(distribution, patch)
    call sparse_collective_root_matches_3d( &
      distribution, root, local_ok, global_ok)
    if (.not. global_ok) return
    call sparse_collective_layout_matches_3d( &
      distribution, patch, hierarchy%nvar, local_ok, global_ok)
    if (.not. global_ok) return
    call gather_rank4_sparse_slabs_3d( &
      distribution%comm, distribution%rank, root, &
      distribution%coarse_counts, distribution%coarse_displacements, &
      hierarchy%coarse_state, global_coarse_state, local_ok)
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    call gather_rank3_sparse_slabs_3d( &
      distribution%comm, distribution%rank, root, &
      distribution%coarse_counts, distribution%coarse_displacements, &
      hierarchy%coarse_temperature, global_coarse_temperature, local_ok)
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    call gather_rank4_sparse_slabs_3d( &
      distribution%comm, distribution%rank, root, &
      distribution%fine_counts, distribution%fine_displacements, &
      hierarchy%fine_state, global_fine_state, local_ok)
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    if (.not. global_ok) return
    call gather_rank3_sparse_slabs_3d( &
      distribution%comm, distribution%rank, root, &
      distribution%fine_counts, distribution%fine_displacements, &
      hierarchy%fine_temperature, global_fine_temperature, local_ok)
    call sparse_collective_logical_and_3d( &
      distribution%comm, local_ok, global_ok)
    ok = global_ok
  end subroutine gather_mpi_amr_sparse_hierarchy_3d

  subroutine scatter_rank4_sparse_slabs_3d( &
      comm, rank, root, slab_counts, slab_displacements, global_values, &
      local_values, ok)
    type(MPI_Comm), intent(in) :: comm
    integer, intent(in) :: rank, root
    integer, intent(in) :: slab_counts(:), slab_displacements(:)
    real(dp), intent(in) :: global_values(:, :, :, :)
    real(dp), intent(out) :: local_values(:, :, :, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: send_values(:), receive_values(:)
    integer, allocatable :: value_counts(:), value_displacements(:)
    integer :: rank_index, plane_values, send_count, ierr
    integer :: first_i, last_i, offset

    plane_values = size(local_values, 1) * size(local_values, 3) * &
      size(local_values, 4)
    allocate(value_counts(size(slab_counts)))
    allocate(value_displacements(size(slab_displacements)))
    value_counts = slab_counts * plane_values
    value_displacements = slab_displacements * plane_values
    send_count = sum(value_counts)
    if (rank == root) then
      allocate(send_values(max(1, send_count)))
      do rank_index = 1, size(slab_counts)
        if (slab_counts(rank_index) == 0) cycle
        first_i = slab_displacements(rank_index) + 1
        last_i = first_i + slab_counts(rank_index) - 1
        offset = value_displacements(rank_index)
        send_values(offset + 1:offset + value_counts(rank_index)) = &
          reshape(global_values(:, first_i:last_i, :, :), &
            [value_counts(rank_index)])
      end do
    else
      allocate(send_values(1))
      send_values = 0.0_dp
    end if
    allocate(receive_values(max(1, size(local_values))))
    receive_values = 0.0_dp
    call MPI_Scatterv( &
      send_values, value_counts, value_displacements, MPI_DOUBLE_PRECISION, &
      receive_values, size(local_values), MPI_DOUBLE_PRECISION, root, &
      comm, ierr)
    ok = ierr == MPI_SUCCESS
    if (.not. ok) return
    if (size(local_values) > 0) then
      local_values = reshape( &
        receive_values(:size(local_values)), shape(local_values))
    end if
  end subroutine scatter_rank4_sparse_slabs_3d

  subroutine scatter_rank3_sparse_slabs_3d( &
      comm, rank, root, slab_counts, slab_displacements, global_values, &
      local_values, ok)
    type(MPI_Comm), intent(in) :: comm
    integer, intent(in) :: rank, root
    integer, intent(in) :: slab_counts(:), slab_displacements(:)
    real(dp), intent(in) :: global_values(:, :, :)
    real(dp), intent(out) :: local_values(:, :, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: send_values(:), receive_values(:)
    integer, allocatable :: value_counts(:), value_displacements(:)
    integer :: rank_index, plane_values, send_count, ierr
    integer :: first_i, last_i, offset

    plane_values = size(local_values, 2) * size(local_values, 3)
    allocate(value_counts(size(slab_counts)))
    allocate(value_displacements(size(slab_displacements)))
    value_counts = slab_counts * plane_values
    value_displacements = slab_displacements * plane_values
    send_count = sum(value_counts)
    if (rank == root) then
      allocate(send_values(max(1, send_count)))
      do rank_index = 1, size(slab_counts)
        if (slab_counts(rank_index) == 0) cycle
        first_i = slab_displacements(rank_index) + 1
        last_i = first_i + slab_counts(rank_index) - 1
        offset = value_displacements(rank_index)
        send_values(offset + 1:offset + value_counts(rank_index)) = &
          reshape(global_values(first_i:last_i, :, :), &
            [value_counts(rank_index)])
      end do
    else
      allocate(send_values(1))
      send_values = 0.0_dp
    end if
    allocate(receive_values(max(1, size(local_values))))
    receive_values = 0.0_dp
    call MPI_Scatterv( &
      send_values, value_counts, value_displacements, MPI_DOUBLE_PRECISION, &
      receive_values, size(local_values), MPI_DOUBLE_PRECISION, root, &
      comm, ierr)
    ok = ierr == MPI_SUCCESS
    if (.not. ok) return
    if (size(local_values) > 0) then
      local_values = reshape( &
        receive_values(:size(local_values)), shape(local_values))
    end if
  end subroutine scatter_rank3_sparse_slabs_3d

  subroutine gather_rank4_sparse_slabs_3d( &
      comm, rank, root, slab_counts, slab_displacements, local_values, &
      global_values, ok)
    type(MPI_Comm), intent(in) :: comm
    integer, intent(in) :: rank, root
    integer, intent(in) :: slab_counts(:), slab_displacements(:)
    real(dp), intent(in) :: local_values(:, :, :, :)
    real(dp), allocatable, intent(out) :: global_values(:, :, :, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: send_values(:), receive_values(:)
    integer, allocatable :: value_counts(:), value_displacements(:)
    integer :: rank_index, plane_values, receive_count, ierr
    integer :: first_i, last_i, offset, global_nx

    plane_values = size(local_values, 1) * size(local_values, 3) * &
      size(local_values, 4)
    allocate(value_counts(size(slab_counts)))
    allocate(value_displacements(size(slab_displacements)))
    value_counts = slab_counts * plane_values
    value_displacements = slab_displacements * plane_values
    receive_count = sum(value_counts)
    allocate(send_values(max(1, size(local_values))))
    send_values = 0.0_dp
    if (size(local_values) > 0) send_values(:size(local_values)) = &
      reshape(local_values, [size(local_values)])
    allocate(receive_values(max(1, receive_count)))
    receive_values = 0.0_dp
    call MPI_Gatherv( &
      send_values, size(local_values), MPI_DOUBLE_PRECISION, receive_values, &
      value_counts, value_displacements, MPI_DOUBLE_PRECISION, root, &
      comm, ierr)
    ok = ierr == MPI_SUCCESS
    if (.not. ok) return
    if (rank == root) then
      global_nx = sum(slab_counts)
      allocate(global_values( &
        size(local_values, 1), global_nx, &
        size(local_values, 3), size(local_values, 4)))
      do rank_index = 1, size(slab_counts)
        if (slab_counts(rank_index) == 0) cycle
        first_i = slab_displacements(rank_index) + 1
        last_i = first_i + slab_counts(rank_index) - 1
        offset = value_displacements(rank_index)
        global_values(:, first_i:last_i, :, :) = reshape( &
          receive_values(offset + 1:offset + value_counts(rank_index)), &
          [size(local_values, 1), slab_counts(rank_index), &
            size(local_values, 3), size(local_values, 4)])
      end do
    else
      allocate(global_values(0, 0, 0, 0))
    end if
  end subroutine gather_rank4_sparse_slabs_3d

  subroutine gather_rank3_sparse_slabs_3d( &
      comm, rank, root, slab_counts, slab_displacements, local_values, &
      global_values, ok)
    type(MPI_Comm), intent(in) :: comm
    integer, intent(in) :: rank, root
    integer, intent(in) :: slab_counts(:), slab_displacements(:)
    real(dp), intent(in) :: local_values(:, :, :)
    real(dp), allocatable, intent(out) :: global_values(:, :, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: send_values(:), receive_values(:)
    integer, allocatable :: value_counts(:), value_displacements(:)
    integer :: rank_index, plane_values, receive_count, ierr
    integer :: first_i, last_i, offset, global_nx

    plane_values = size(local_values, 2) * size(local_values, 3)
    allocate(value_counts(size(slab_counts)))
    allocate(value_displacements(size(slab_displacements)))
    value_counts = slab_counts * plane_values
    value_displacements = slab_displacements * plane_values
    receive_count = sum(value_counts)
    allocate(send_values(max(1, size(local_values))))
    send_values = 0.0_dp
    if (size(local_values) > 0) send_values(:size(local_values)) = &
      reshape(local_values, [size(local_values)])
    allocate(receive_values(max(1, receive_count)))
    receive_values = 0.0_dp
    call MPI_Gatherv( &
      send_values, size(local_values), MPI_DOUBLE_PRECISION, receive_values, &
      value_counts, value_displacements, MPI_DOUBLE_PRECISION, root, &
      comm, ierr)
    ok = ierr == MPI_SUCCESS
    if (.not. ok) return
    if (rank == root) then
      global_nx = sum(slab_counts)
      allocate(global_values( &
        global_nx, size(local_values, 2), size(local_values, 3)))
      do rank_index = 1, size(slab_counts)
        if (slab_counts(rank_index) == 0) cycle
        first_i = slab_displacements(rank_index) + 1
        last_i = first_i + slab_counts(rank_index) - 1
        offset = value_displacements(rank_index)
        global_values(first_i:last_i, :, :) = reshape( &
          receive_values(offset + 1:offset + value_counts(rank_index)), &
          [slab_counts(rank_index), &
            size(local_values, 2), size(local_values, 3)])
      end do
    else
      allocate(global_values(0, 0, 0))
    end if
  end subroutine gather_rank3_sparse_slabs_3d

  subroutine sparse_collective_contract_matches_3d( &
      distribution, species, patch, local_reals, local_text, &
      local_valid, ok)
    type(mpi_amr_sparse_distribution_3d), intent(in) :: distribution
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: local_reals(:)
    character(len=*), intent(in) :: local_text
    logical, intent(in) :: local_valid
    logical, intent(out) :: ok

    real(dp), allocatable :: root_reals(:)
    real(dp) :: local_species_reals( &
      sparse_nasa7_real_count_3d, sparse_maximum_species_count_3d)
    real(dp) :: root_species_reals( &
      sparse_nasa7_real_count_3d, sparse_maximum_species_count_3d)
    integer :: local_patch(11), root_patch(11), species_index, ierr
    character(len=24) :: &
      local_species_names(sparse_maximum_species_count_3d)
    character(len=24) :: &
      root_species_names(sparse_maximum_species_count_3d)
    character(len=256) :: root_text
    logical :: matches, global_matches, layout_matches, collective_ok

    call sparse_collective_layout_matches_3d( &
      distribution, patch, reactive_nvar(size(species)), local_valid, &
      layout_matches)
    if (.not. layout_matches) then
      ok = .false.
      return
    end if

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
    if (size(species) <= sparse_maximum_species_count_3d) then
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
    root_text = ""
    if (distribution%rank == 0) root_text = trim(local_text)
    call MPI_Bcast(root_reals, size(root_reals), MPI_DOUBLE_PRECISION, 0, &
      distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, collective_ok)
    if (.not. collective_ok) return
    call MPI_Bcast(root_patch, size(root_patch), MPI_INTEGER, 0, &
      distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, collective_ok)
    if (.not. collective_ok) return
    call MPI_Bcast(root_species_reals, size(root_species_reals), &
      MPI_DOUBLE_PRECISION, 0, distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, collective_ok)
    if (.not. collective_ok) return
    call MPI_Bcast(root_species_names, &
      len(root_species_names) * size(root_species_names), MPI_CHARACTER, &
      0, distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, collective_ok)
    if (.not. collective_ok) return
    call MPI_Bcast(root_text, len(root_text), MPI_CHARACTER, 0, &
      distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, collective_ok)
    if (.not. collective_ok) return
    matches = local_valid .and. &
      size(species) <= sparse_maximum_species_count_3d .and. &
      all(sparse_identical_real_bits_3d(local_reals, root_reals)) .and. &
      all(local_patch == root_patch) .and. &
      all(sparse_identical_real_bits_3d( &
        local_species_reals, root_species_reals)) .and. &
      all(local_species_names == root_species_names) .and. &
      trim(local_text) == trim(root_text)
    call sparse_collective_logical_and_3d( &
      distribution%comm, matches, global_matches)
    ok = global_matches
  end subroutine sparse_collective_contract_matches_3d

  subroutine sparse_collective_transport_matches_3d( &
      distribution, transport, local_valid, ok)
    type(mpi_amr_sparse_distribution_3d), intent(in) :: distribution
    type(gas_transport_species), intent(in) :: transport(:)
    logical, intent(in) :: local_valid
    logical, intent(out) :: ok

    real(dp), allocatable :: local_reals(:, :), root_reals(:, :)
    integer, allocatable :: local_geometry(:), root_geometry(:)
    character(len=24), allocatable :: local_names(:), root_names(:)
    integer :: local_count, minimum_count, maximum_count
    integer :: species_index, ierr
    logical :: counts_match, local_ok, matches, collective_ok

    ok = .false.
    local_count = size(transport)
    call MPI_Allreduce( &
      local_count, minimum_count, 1, MPI_INTEGER, MPI_MIN, &
      distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, collective_ok)
    if (.not. collective_ok) return
    call MPI_Allreduce( &
      local_count, maximum_count, 1, MPI_INTEGER, MPI_MAX, &
      distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, collective_ok)
    if (.not. collective_ok) return
    counts_match = local_valid .and. minimum_count >= 1 .and. &
      minimum_count == maximum_count
    call sparse_collective_logical_and_3d( &
      distribution%comm, counts_match, collective_ok)
    if (.not. collective_ok) return

    allocate(local_reals(5, local_count), root_reals(5, local_count))
    allocate(local_geometry(local_count), root_geometry(local_count))
    allocate(local_names(local_count), root_names(local_count))
    local_ok = local_valid
    do species_index = 1, local_count
      local_names(species_index) = transport(species_index)%name
      local_geometry(species_index) = transport(species_index)%geometry
      local_reals(:, species_index) = [ &
        transport(species_index)%well_depth, &
        transport(species_index)%diameter, &
        transport(species_index)%dipole, &
        transport(species_index)%polarizability, &
        transport(species_index)%rotational_relaxation]
    end do
    if (local_ok) local_ok = all(ieee_is_finite(local_reals))
    if (local_ok) then
      do species_index = 1, local_count
        local_ok = len_trim(local_names(species_index)) > 0 .and. &
          local_geometry(species_index) >= 0 .and. &
          local_geometry(species_index) <= 2 .and. &
          local_reals(1, species_index) > 0.0_dp .and. &
          local_reals(2, species_index) > 0.0_dp .and. &
          all(local_reals(3:5, species_index) >= 0.0_dp)
        if (.not. local_ok) exit
      end do
    end if
    root_reals = local_reals
    root_geometry = local_geometry
    root_names = local_names
    call MPI_Bcast( &
      root_reals, size(root_reals), MPI_DOUBLE_PRECISION, 0, &
      distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, collective_ok)
    if (.not. collective_ok) return
    call MPI_Bcast( &
      root_geometry, size(root_geometry), MPI_INTEGER, 0, &
      distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, collective_ok)
    if (.not. collective_ok) return
    call MPI_Bcast( &
      root_names, len(root_names) * size(root_names), MPI_CHARACTER, 0, &
      distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, collective_ok)
    if (.not. collective_ok) return
    matches = local_ok .and. &
      all(sparse_identical_real_bits_3d(local_reals, root_reals)) .and. &
      all(local_geometry == root_geometry) .and. all(local_names == root_names)
    call sparse_collective_logical_and_3d( &
      distribution%comm, matches, ok)
  end subroutine sparse_collective_transport_matches_3d

  subroutine sparse_collective_reactions_match_3d( &
      distribution, reactions, nspecies, local_valid, ok)
    type(mpi_amr_sparse_distribution_3d), intent(in) :: distribution
    type(elementary_reaction), intent(in) :: reactions(:)
    integer, intent(in) :: nspecies
    logical, intent(in) :: local_valid
    logical, intent(out) :: ok

    real(dp), allocatable :: local_reals(:, :), root_reals(:, :)
    integer, allocatable :: local_integers(:, :), root_integers(:, :)
    character(len=128), allocatable :: local_equations(:), root_equations(:)
    integer :: local_count, minimum_count, maximum_count
    integer :: reaction_index, real_index, real_count, ierr
    logical :: local_ok, reaction_ok, counts_match, matches, collective_ok

    ok = .false.
    local_count = size(reactions)
    call MPI_Allreduce( &
      local_count, minimum_count, 1, MPI_INTEGER, MPI_MIN, &
      distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, collective_ok)
    if (.not. collective_ok) return
    call MPI_Allreduce( &
      local_count, maximum_count, 1, MPI_INTEGER, MPI_MAX, &
      distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, collective_ok)
    if (.not. collective_ok) return
    counts_match = local_valid .and. nspecies >= 1 .and. &
      minimum_count >= 1 .and. minimum_count == maximum_count
    call sparse_collective_logical_and_3d( &
      distribution%comm, counts_match, collective_ok)
    if (.not. collective_ok) return

    real_count = 3 * nspecies + 13
    allocate(local_reals(real_count, local_count))
    allocate(root_reals(real_count, local_count))
    allocate(local_integers(4, local_count))
    allocate(root_integers(4, local_count))
    allocate(local_equations(local_count), root_equations(local_count))
    local_reals = 0.0_dp
    local_integers = 0
    local_equations = ""
    local_ok = local_valid
    do reaction_index = 1, local_count
      reaction_ok = valid_elementary_reaction( &
        reactions(reaction_index), nspecies)
      local_ok = local_ok .and. reaction_ok
      local_integers(:, reaction_index) = [ &
        reactions(reaction_index)%kind, &
        merge(1, 0, reactions(reaction_index)%reversible), &
        merge(1, 0, reactions(reaction_index)%troe%enabled), &
        merge(1, 0, allocated( &
          reactions(reaction_index)%third_body_efficiencies))]
      local_equations(reaction_index) = reactions(reaction_index)%equation
      real_index = 1
      if (allocated(reactions(reaction_index)%reactant_stoich)) then
        if (size(reactions(reaction_index)%reactant_stoich) == nspecies) &
          local_reals(real_index:real_index + nspecies - 1, &
            reaction_index) = reactions(reaction_index)%reactant_stoich
      end if
      real_index = real_index + nspecies
      if (allocated(reactions(reaction_index)%product_stoich)) then
        if (size(reactions(reaction_index)%product_stoich) == nspecies) &
          local_reals(real_index:real_index + nspecies - 1, &
            reaction_index) = reactions(reaction_index)%product_stoich
      end if
      real_index = real_index + nspecies
      local_reals(real_index:real_index + 8, reaction_index) = [ &
        reactions(reaction_index)%forward_rate%pre_exponential, &
        reactions(reaction_index)%forward_rate%temperature_exponent, &
        reactions(reaction_index)%forward_rate%activation_energy, &
        reactions(reaction_index)%low_pressure_rate%pre_exponential, &
        reactions(reaction_index)%low_pressure_rate%temperature_exponent, &
        reactions(reaction_index)%low_pressure_rate%activation_energy, &
        reactions(reaction_index)%high_pressure_rate%pre_exponential, &
        reactions(reaction_index)%high_pressure_rate%temperature_exponent, &
        reactions(reaction_index)%high_pressure_rate%activation_energy]
      real_index = real_index + 9
      if (allocated( &
          reactions(reaction_index)%third_body_efficiencies)) then
        if (size(reactions(reaction_index)%third_body_efficiencies) == &
            nspecies) local_reals( &
          real_index:real_index + nspecies - 1, reaction_index) = &
            reactions(reaction_index)%third_body_efficiencies
      end if
      real_index = real_index + nspecies
      local_reals(real_index:real_index + 3, reaction_index) = [ &
        reactions(reaction_index)%troe%alpha, &
        reactions(reaction_index)%troe%temperature_3, &
        reactions(reaction_index)%troe%temperature_1, &
        reactions(reaction_index)%troe%temperature_2]
    end do
    root_reals = local_reals
    root_integers = local_integers
    root_equations = local_equations
    call MPI_Bcast( &
      root_reals, size(root_reals), MPI_DOUBLE_PRECISION, 0, &
      distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, collective_ok)
    if (.not. collective_ok) return
    call MPI_Bcast( &
      root_integers, size(root_integers), MPI_INTEGER, 0, &
      distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, collective_ok)
    if (.not. collective_ok) return
    call MPI_Bcast( &
      root_equations, len(root_equations) * size(root_equations), &
      MPI_CHARACTER, 0, distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, collective_ok)
    if (.not. collective_ok) return
    matches = local_ok .and. &
      all(sparse_identical_real_bits_3d(local_reals, root_reals)) .and. &
      all(local_integers == root_integers) .and. &
      all(local_equations == root_equations)
    call sparse_collective_logical_and_3d( &
      distribution%comm, matches, ok)
  end subroutine sparse_collective_reactions_match_3d

  subroutine sparse_collective_root_matches_3d( &
      distribution, local_root, local_valid, ok)
    type(mpi_amr_sparse_distribution_3d), intent(in) :: distribution
    integer, intent(in) :: local_root
    logical, intent(in) :: local_valid
    logical, intent(out) :: ok

    integer :: minimum_root, maximum_root, ierr
    logical :: roots_match, collective_ok

    call MPI_Allreduce(local_root, minimum_root, 1, MPI_INTEGER, MPI_MIN, &
      distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, collective_ok)
    if (.not. collective_ok) return
    call MPI_Allreduce(local_root, maximum_root, 1, MPI_INTEGER, MPI_MAX, &
      distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, collective_ok)
    if (.not. collective_ok) return
    roots_match = local_valid .and. minimum_root == maximum_root .and. &
      minimum_root >= 0 .and. minimum_root < distribution%nranks
    call sparse_collective_logical_and_3d( &
      distribution%comm, roots_match, ok)
  end subroutine sparse_collective_root_matches_3d

  subroutine sparse_collective_layout_matches_3d( &
      distribution, patch, nvar, local_valid, ok)
    type(mpi_amr_sparse_distribution_3d), intent(in) :: distribution
    type(amr_patch_3d), intent(in) :: patch
    integer, intent(in) :: nvar
    logical, intent(in) :: local_valid
    logical, intent(out) :: ok

    integer, allocatable :: local_layout(:), root_layout(:)
    integer :: comm_rank, comm_size, first, last, ierr
    logical :: arrays_valid, matches, collective_ok

    ok = .false.
    call MPI_Comm_rank(distribution%comm, comm_rank, ierr)
    if (ierr /= MPI_SUCCESS) return
    call MPI_Comm_size(distribution%comm, comm_size, ierr)
    if (ierr /= MPI_SUCCESS) return
    allocate(local_layout(12 + 4 * comm_size), root_layout(12 + 4 * comm_size))
    local_layout = -huge(1)
    local_layout(1:12) = [ &
      patch%coarse_nx, patch%coarse_ny, patch%coarse_nz, &
      patch%coarse_i_lower, patch%coarse_i_upper, &
      patch%coarse_j_lower, patch%coarse_j_upper, &
      patch%coarse_k_lower, patch%coarse_k_upper, patch%refinement_ratio, &
      distribution%nranks, nvar]
    arrays_valid = allocated(distribution%coarse_counts) .and. &
      allocated(distribution%coarse_displacements) .and. &
      allocated(distribution%fine_counts) .and. &
      allocated(distribution%fine_displacements)
    if (arrays_valid) then
      arrays_valid = size(distribution%coarse_counts) == comm_size .and. &
        size(distribution%coarse_displacements) == comm_size .and. &
        size(distribution%fine_counts) == comm_size .and. &
        size(distribution%fine_displacements) == comm_size
    end if
    if (arrays_valid) then
      first = 13
      last = first + comm_size - 1
      local_layout(first:last) = distribution%coarse_counts
      first = last + 1
      last = first + comm_size - 1
      local_layout(first:last) = distribution%coarse_displacements
      first = last + 1
      last = first + comm_size - 1
      local_layout(first:last) = distribution%fine_counts
      first = last + 1
      last = first + comm_size - 1
      local_layout(first:last) = distribution%fine_displacements
    end if
    root_layout = local_layout
    call MPI_Bcast( &
      root_layout, size(root_layout), MPI_INTEGER, 0, &
      distribution%comm, ierr)
    call sparse_collective_mpi_success_3d( &
      distribution%comm, ierr, collective_ok)
    if (.not. collective_ok) return
    matches = local_valid .and. arrays_valid .and. &
      distribution%rank == comm_rank .and. &
      distribution%nranks == comm_size .and. &
      all(local_layout == root_layout)
    call sparse_collective_logical_and_3d(distribution%comm, matches, ok)
  end subroutine sparse_collective_layout_matches_3d

  pure elemental logical function sparse_identical_real_bits_3d( &
      left, right) result(identical)
    real(dp), intent(in) :: left, right

    identical = transfer(left, 0_int64) == transfer(right, 0_int64)
  end function sparse_identical_real_bits_3d

  subroutine sparse_collective_logical_and_3d(comm, local_value, global_value)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(in) :: local_value
    logical, intent(out) :: global_value

    integer :: ierr

    call MPI_Allreduce( &
      local_value, global_value, 1, MPI_LOGICAL, MPI_LAND, comm, ierr)
    if (ierr /= MPI_SUCCESS) global_value = .false.
  end subroutine sparse_collective_logical_and_3d

  subroutine sparse_collective_mpi_success_3d( &
      comm, local_ierr, global_success)
    type(MPI_Comm), intent(in) :: comm
    integer, intent(in) :: local_ierr
    logical, intent(out) :: global_success

    logical :: local_success
    integer :: ierr

    local_success = local_ierr == MPI_SUCCESS
    call MPI_Allreduce( &
      local_success, global_success, 1, MPI_LOGICAL, MPI_LAND, comm, ierr)
    if (ierr /= MPI_SUCCESS) global_success = .false.
  end subroutine sparse_collective_mpi_success_3d

end module mpi_amr_sparse_reactive_3d_mod
