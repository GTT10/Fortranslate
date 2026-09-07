program test_amr_reactive_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, &
    ieee_quiet_nan
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mesh_3d_mod, only: uniform_cell_centers_3d
  use reactive_1d_mod, only: reactive_nvar, reactive_nprim
  use simulation_config_reactive_3d_mod, only: reactive_3d_config
  use reactive_entropy_wave_3d_problem_mod, only: &
    initialize_reactive_problem_3d, reactive_entropy_wave_density_3d
  use reactive_3d_mod, only: &
    recover_reactive_temperatures_3d, reactive_extrema_3d
  use amr_hierarchy_3d_mod, only: &
    amr_patch_3d, initialize_amr_patch_3d, restrict_average_3d, &
    average_down_3d, composite_integrals_amr_3d
  use amr_reactive_3d_mod, only: &
    compute_amr_reactive_cfl_timestep_3d, &
    advance_amr_reactive_hydro_3d, accumulate_fine_interface_fluxes_3d, &
    reflux_coarse_3d, interpolate_coarse_fine_ghost_plm_3d
  implicit none

  integer, parameter :: n = 8, ratio = 2
  real(dp), parameter :: final_time = 1.0e-6_dp
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_3d) :: patch
  real(dp) :: conservation_error, synchronization_error
  real(dp) :: maximum_reflux, closure_error
  real(dp) :: density_error, plm_density_error
  real(dp) :: plm_conservation_error, plm_synchronization_error
  real(dp) :: plm_maximum_reflux, plm_closure_error
  logical :: ok

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "elementary thermodynamics load")
  call initialize_amr_patch_3d( &
    n, n, n, 3, 6, 3, 6, 3, 6, ratio, patch, ok)
  call require(ok .and. patch%is_strictly_interior(), &
    "strictly interior two-level patch")
  call run_entropy_wave( &
    species, patch, "pcm", conservation_error, synchronization_error, &
    maximum_reflux, closure_error, density_error, ok)
  call require(ok, "nonuniform AMR advance")
  call require(conservation_error <= 2.0e-11_dp, &
    "refluxed composite conservation")
  call require(synchronization_error <= 5.0e-13_dp, &
    "covered coarse average-down synchronization")
  call require(maximum_reflux > 1.0e-12_dp, &
    "nonuniform flow exercises reflux")
  call require(closure_error <= 5.0e-11_dp, &
    "coarse and fine species closure")
  call run_entropy_wave( &
    species, patch, "characteristic_plm", plm_conservation_error, &
    plm_synchronization_error, plm_maximum_reflux, plm_closure_error, &
    plm_density_error, ok)
  call require(ok, "characteristic PLM nonuniform AMR advance")
  call require(plm_conservation_error <= 2.0e-11_dp, &
    "characteristic PLM refluxed composite conservation")
  call require(plm_synchronization_error <= 5.0e-13_dp, &
    "characteristic PLM average-down synchronization")
  call require(plm_maximum_reflux > 1.0e-12_dp, &
    "characteristic PLM exercises reflux")
  call require(plm_closure_error <= 5.0e-11_dp, &
    "characteristic PLM species closure")
  call require(plm_density_error <= 0.70_dp * density_error, &
    "characteristic PLM fine-grid accuracy improvement")
  call run_uniform_invariance(species, patch, "pcm", "mc", ok)
  call require(ok, "uniform PCM AMR invariance")
  call run_uniform_invariance( &
    species, patch, "characteristic_plm", "mc", ok)
  call require(ok, "uniform characteristic PLM AMR invariance")
  call run_uniform_invariance( &
    species, patch, "characteristic_plm", "minmod", ok)
  call require(ok, "uniform minmod characteristic PLM AMR invariance")
  call run_public_api_validation_tests(species, patch, ok)
  call require(ok, "public AMR API validation and rollback")
  write(*, '(a,es24.16)') "composite conservation error=", conservation_error
  write(*, '(a,es24.16)') "average-down synchronization error=", &
    synchronization_error
  write(*, '(a,es24.16)') "maximum reflux correction=", maximum_reflux
  write(*, '(a,es24.16)') "maximum species closure error=", closure_error
  write(*, '(a,es24.16)') "PCM fine density L1 error=", density_error
  write(*, '(a,es24.16)') &
    "PLM composite conservation error=", plm_conservation_error
  write(*, '(a,es24.16)') &
    "PLM average-down synchronization error=", plm_synchronization_error
  write(*, '(a,es24.16)') &
    "PLM maximum reflux correction=", plm_maximum_reflux
  write(*, '(a,es24.16)') &
    "PLM maximum species closure error=", plm_closure_error
  write(*, '(a,es24.16)') &
    "PLM fine density L1 error=", plm_density_error
  write(*, '(a)') "test_amr_reactive_3d: PASS"

contains

  subroutine run_public_api_validation_tests(species, patch, ok_out)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    logical, intent(out) :: ok_out

    type(amr_patch_3d) :: invalid_patch
    real(dp), allocatable :: coarse_state(:, :, :, :)
    real(dp), allocatable :: coarse_temperature(:, :, :)
    real(dp), allocatable :: fine_state(:, :, :, :)
    real(dp), allocatable :: fine_temperature(:, :, :)
    real(dp), allocatable :: coarse_start(:, :, :, :)
    real(dp), allocatable :: coarse_end(:, :, :, :)
    real(dp), allocatable :: coarse_start_temperature(:, :, :)
    real(dp), allocatable :: coarse_end_temperature(:, :, :)
    real(dp), allocatable :: bad_coarse_start(:, :, :, :)
    real(dp), allocatable :: primitive(:), ghost_state(:)
    real(dp), allocatable :: face_flux_x(:, :, :, :)
    real(dp), allocatable :: face_flux_y(:, :, :, :)
    real(dp), allocatable :: face_flux_z(:, :, :, :)
    real(dp), allocatable :: x_lower(:, :, :), x_upper(:, :, :)
    real(dp), allocatable :: y_lower(:, :, :), y_upper(:, :, :)
    real(dp), allocatable :: z_lower(:, :, :), z_upper(:, :, :)
    real(dp), allocatable :: saved_x_lower(:, :, :), saved_x_upper(:, :, :)
    real(dp), allocatable :: saved_y_lower(:, :, :), saved_y_upper(:, :, :)
    real(dp), allocatable :: saved_z_lower(:, :, :), saved_z_upper(:, :, :)
    real(dp), allocatable :: bad_x_lower(:, :, :), saved_bad_x_lower(:, :, :)
    real(dp), allocatable :: coarse_flux_x(:, :, :, :)
    real(dp), allocatable :: coarse_flux_y(:, :, :, :)
    real(dp), allocatable :: coarse_flux_z(:, :, :, :)
    real(dp), allocatable :: fine_x_lower(:, :, :), fine_x_upper(:, :, :)
    real(dp), allocatable :: fine_y_lower(:, :, :), fine_y_upper(:, :, :)
    real(dp), allocatable :: fine_z_lower(:, :, :), fine_z_upper(:, :, :)
    real(dp), allocatable :: saved_coarse_state(:, :, :, :)
    real(dp), allocatable :: saved_bad_coarse_state(:, :, :, :)
    real(dp) :: ghost_temperature, maximum_correction
    integer :: nvar, nprim, nx_coarse, ny_coarse, nz_coarse
    integer :: nx_fine, ny_fine, nz_fine
    integer :: covered_nx, covered_ny, covered_nz
    logical :: local_ok

    ok_out = .false.
    call initialize_public_api_hierarchy( &
      species, patch, coarse_state, coarse_temperature, fine_state, &
      fine_temperature, local_ok)
    if (.not. local_ok) return
    nvar = reactive_nvar(size(species))
    nprim = reactive_nprim(size(species))
    nx_coarse = size(coarse_state, 2)
    ny_coarse = size(coarse_state, 3)
    nz_coarse = size(coarse_state, 4)
    nx_fine = patch%fine_nx()
    ny_fine = patch%fine_ny()
    nz_fine = patch%fine_nz()
    covered_nx = patch%coarse_i_upper - patch%coarse_i_lower + 1
    covered_ny = patch%coarse_j_upper - patch%coarse_j_lower + 1
    covered_nz = patch%coarse_k_upper - patch%coarse_k_lower + 1
    allocate(coarse_start, source=coarse_state)
    allocate(coarse_end, source=coarse_state)
    allocate(coarse_start_temperature, source=coarse_temperature)
    allocate(coarse_end_temperature, source=coarse_temperature)
    allocate(primitive(nprim), ghost_state(nvar))

    call interpolate_coarse_fine_ghost_plm_3d( &
      species, coarse_start, coarse_start_temperature, coarse_end, &
      coarse_end_temperature, patch%coarse_i_lower, patch%coarse_j_lower, &
      patch%coarse_k_lower, 0, 1, 1, patch%refinement_ratio, 0.5_dp, "mc", &
      primitive, ghost_state, ghost_temperature, local_ok)
    call require(local_ok, "valid public PLM ghost call")
    call interpolate_coarse_fine_ghost_plm_3d( &
      species, coarse_start, coarse_start_temperature, coarse_end, &
      coarse_end_temperature, patch%coarse_i_lower, patch%coarse_j_lower, &
      patch%coarse_k_lower, 0, 1, 1, 0, 0.5_dp, "mc", primitive, &
      ghost_state, ghost_temperature, local_ok)
    call require(.not. local_ok, "public PLM ghost ratio-zero rejection")
    call interpolate_coarse_fine_ghost_plm_3d( &
      species, coarse_start, coarse_start_temperature, coarse_end, &
      coarse_end_temperature, 0, patch%coarse_j_lower, patch%coarse_k_lower, &
      0, 1, 1, patch%refinement_ratio, 0.5_dp, "mc", primitive, ghost_state, &
      ghost_temperature, local_ok)
    call require(.not. local_ok, "public PLM ghost bad-index rejection")
    call interpolate_coarse_fine_ghost_plm_3d( &
      species, coarse_start, coarse_start_temperature, coarse_end, &
      coarse_end_temperature, patch%coarse_i_lower, patch%coarse_j_lower, &
      patch%coarse_k_lower, -2, 1, 1, patch%refinement_ratio, 0.5_dp, "mc", &
      primitive, ghost_state, ghost_temperature, local_ok)
    call require(.not. local_ok, "public PLM ghost local-index rejection")
    call interpolate_coarse_fine_ghost_plm_3d( &
      species, coarse_start, coarse_start_temperature, coarse_end, &
      coarse_end_temperature, patch%coarse_i_lower, patch%coarse_j_lower, &
      patch%coarse_k_lower, 0, 1, 1, patch%refinement_ratio, 0.5_dp, &
      "invalid", primitive, ghost_state, ghost_temperature, local_ok)
    call require(.not. local_ok, "public PLM ghost limiter rejection")
    call interpolate_coarse_fine_ghost_plm_3d( &
      species, coarse_start, coarse_start_temperature, coarse_end, &
      coarse_end_temperature, patch%coarse_i_lower, patch%coarse_j_lower, &
      patch%coarse_k_lower, 0, 1, 1, patch%refinement_ratio, -0.5_dp, "mc", &
      primitive, ghost_state, ghost_temperature, local_ok)
    call require(.not. local_ok, "public PLM ghost alpha rejection")
    allocate(bad_coarse_start(nvar, nx_coarse - 1, ny_coarse, nz_coarse))
    bad_coarse_start = coarse_start(:, 1:nx_coarse - 1, :, :)
    call interpolate_coarse_fine_ghost_plm_3d( &
      species, bad_coarse_start, coarse_start_temperature, coarse_end, &
      coarse_end_temperature, patch%coarse_i_lower, patch%coarse_j_lower, &
      patch%coarse_k_lower, 0, 1, 1, patch%refinement_ratio, 0.5_dp, "mc", &
      primitive, ghost_state, ghost_temperature, local_ok)
    call require(.not. local_ok, "public PLM ghost shape rejection")

    allocate(face_flux_x(nvar, 0:nx_fine, ny_fine, nz_fine))
    allocate(face_flux_y(nvar, nx_fine, 0:ny_fine, nz_fine))
    allocate(face_flux_z(nvar, nx_fine, ny_fine, 0:nz_fine))
    allocate(x_lower(nvar, covered_ny, covered_nz), &
      x_upper(nvar, covered_ny, covered_nz))
    allocate(y_lower(nvar, covered_nx, covered_nz), &
      y_upper(nvar, covered_nx, covered_nz))
    allocate(z_lower(nvar, covered_nx, covered_ny), &
      z_upper(nvar, covered_nx, covered_ny))
    face_flux_x = 0.125_dp
    face_flux_y = 0.250_dp
    face_flux_z = 0.375_dp
    x_lower = 1.25_dp
    x_upper = 1.50_dp
    y_lower = 1.75_dp
    y_upper = 2.00_dp
    z_lower = 2.25_dp
    z_upper = 2.50_dp
    call accumulate_fine_interface_fluxes_3d( &
      patch, face_flux_x, face_flux_y, face_flux_z, x_lower, x_upper, &
      y_lower, y_upper, z_lower, z_upper, local_ok)
    call require(local_ok, "valid interface accumulation")
    allocate(saved_x_lower, source=x_lower)
    allocate(saved_x_upper, source=x_upper)
    allocate(saved_y_lower, source=y_lower)
    allocate(saved_y_upper, source=y_upper)
    allocate(saved_z_lower, source=z_lower)
    allocate(saved_z_upper, source=z_upper)
    face_flux_x(1, 0, 1, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
    call accumulate_fine_interface_fluxes_3d( &
      patch, face_flux_x, face_flux_y, face_flux_z, x_lower, x_upper, &
      y_lower, y_upper, z_lower, z_upper, local_ok)
    call require(.not. local_ok .and. all(x_lower == saved_x_lower) .and. &
      all(x_upper == saved_x_upper) .and. all(y_lower == saved_y_lower) .and. &
      all(y_upper == saved_y_upper) .and. all(z_lower == saved_z_lower) .and. &
      all(z_upper == saved_z_upper), &
      "NaN interface flux preserves caller accumulators")
    face_flux_x(1, 0, 1, 1) = 0.125_dp
    allocate(bad_x_lower(nvar, covered_ny - 1, covered_nz))
    allocate(saved_bad_x_lower, source=bad_x_lower)
    bad_x_lower = 3.25_dp
    saved_bad_x_lower = bad_x_lower
    call accumulate_fine_interface_fluxes_3d( &
      patch, face_flux_x, face_flux_y, face_flux_z, bad_x_lower, x_upper, &
      y_lower, y_upper, z_lower, z_upper, local_ok)
    call require(.not. local_ok .and. all(bad_x_lower == saved_bad_x_lower) .and. &
      all(x_upper == saved_x_upper) .and. all(y_lower == saved_y_lower) .and. &
      all(y_upper == saved_y_upper) .and. all(z_lower == saved_z_lower) .and. &
      all(z_upper == saved_z_upper), &
      "bad interface shape preserves caller accumulators")
    invalid_patch = patch
    invalid_patch%refinement_ratio = 0
    call accumulate_fine_interface_fluxes_3d( &
      invalid_patch, face_flux_x, face_flux_y, face_flux_z, x_lower, x_upper, &
      y_lower, y_upper, z_lower, z_upper, local_ok)
    call require(.not. local_ok .and. all(x_lower == saved_x_lower) .and. &
      all(x_upper == saved_x_upper) .and. all(y_lower == saved_y_lower) .and. &
      all(y_upper == saved_y_upper) .and. all(z_lower == saved_z_lower) .and. &
      all(z_upper == saved_z_upper), "zero-ratio interface rejection")

    allocate(coarse_flux_x(nvar, nx_coarse, ny_coarse, nz_coarse))
    allocate(coarse_flux_y(nvar, nx_coarse, ny_coarse, nz_coarse))
    allocate(coarse_flux_z(nvar, nx_coarse, ny_coarse, nz_coarse))
    allocate(fine_x_lower(nvar, covered_ny, covered_nz), &
      fine_x_upper(nvar, covered_ny, covered_nz))
    allocate(fine_y_lower(nvar, covered_nx, covered_nz), &
      fine_y_upper(nvar, covered_nx, covered_nz))
    allocate(fine_z_lower(nvar, covered_nx, covered_ny), &
      fine_z_upper(nvar, covered_nx, covered_ny))
    coarse_flux_x = 0.125_dp
    coarse_flux_y = 0.250_dp
    coarse_flux_z = 0.375_dp
    fine_x_lower = 0.500_dp
    fine_x_upper = 0.625_dp
    fine_y_lower = 0.750_dp
    fine_y_upper = 0.875_dp
    fine_z_lower = 1.000_dp
    fine_z_upper = 1.125_dp
    allocate(saved_coarse_state, source=coarse_state)
    call reflux_coarse_3d( &
      patch, coarse_state, coarse_flux_x, coarse_flux_y, coarse_flux_z, &
      fine_x_lower, fine_x_upper, fine_y_lower, fine_y_upper, fine_z_lower, &
      fine_z_upper, 1.0_dp, 1.0_dp, 1.0_dp, 1.0e-3_dp, &
      maximum_correction, local_ok)
    call require(local_ok .and. maximum_correction > 0.0_dp, &
      "valid coarse reflux")
    coarse_state = saved_coarse_state
    allocate(saved_bad_coarse_state, source=coarse_state)
    coarse_flux_x(1, patch%coarse_i_lower - 1, patch%coarse_j_lower, &
      patch%coarse_k_lower) = ieee_value(0.0_dp, ieee_quiet_nan)
    call reflux_coarse_3d( &
      patch, coarse_state, coarse_flux_x, coarse_flux_y, coarse_flux_z, &
      fine_x_lower, fine_x_upper, fine_y_lower, fine_y_upper, fine_z_lower, &
      fine_z_upper, 1.0_dp, 1.0_dp, 1.0_dp, 1.0e-3_dp, &
      maximum_correction, local_ok)
    call require(.not. local_ok .and. maximum_correction == 0.0_dp .and. &
      all(coarse_state == saved_bad_coarse_state), &
      "NaN coarse reflux preserves caller state")
    coarse_flux_x(1, patch%coarse_i_lower - 1, patch%coarse_j_lower, &
      patch%coarse_k_lower) = 0.125_dp
    call reflux_coarse_3d( &
      patch, coarse_state, coarse_flux_x, coarse_flux_y, coarse_flux_z, &
      bad_x_lower, fine_x_upper, fine_y_lower, fine_y_upper, fine_z_lower, &
      fine_z_upper, 1.0_dp, 1.0_dp, 1.0_dp, 1.0e-3_dp, &
      maximum_correction, local_ok)
    call require(.not. local_ok .and. maximum_correction == 0.0_dp .and. &
      all(coarse_state == saved_bad_coarse_state), &
      "bad fine reflux shape preserves caller state")
    invalid_patch = patch
    invalid_patch%refinement_ratio = 0
    call reflux_coarse_3d( &
      invalid_patch, coarse_state, coarse_flux_x, coarse_flux_y, coarse_flux_z, &
      fine_x_lower, fine_x_upper, fine_y_lower, fine_y_upper, fine_z_lower, &
      fine_z_upper, 1.0_dp, 1.0_dp, 1.0_dp, 1.0e-3_dp, &
      maximum_correction, local_ok)
    call require(.not. local_ok .and. maximum_correction == 0.0_dp .and. &
      all(coarse_state == saved_bad_coarse_state), "zero-ratio reflux rejection")
    ok_out = .true.
  end subroutine run_public_api_validation_tests

  subroutine initialize_public_api_hierarchy( &
      species, patch, coarse_state, coarse_temperature, fine_state, &
      fine_temperature, ok_out)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), allocatable, intent(out) :: coarse_state(:, :, :, :)
    real(dp), allocatable, intent(out) :: coarse_temperature(:, :, :)
    real(dp), allocatable, intent(out) :: fine_state(:, :, :, :)
    real(dp), allocatable, intent(out) :: fine_temperature(:, :, :)
    logical, intent(out) :: ok_out

    type(reactive_3d_config) :: config, fine_config
    real(dp), allocatable :: x(:), y(:), z(:), xf(:), yf(:), zf(:)
    real(dp), allocatable :: mass_fractions(:)
    real(dp) :: dx, dy, dz, density, fine_density
    logical :: local_ok
    integer :: nvar

    ok_out = .false.
    config = reactive_3d_config()
    config%nx = n
    config%ny = n
    config%nz = n
    config%thermo_model = "elementary"
    nvar = reactive_nvar(size(species))
    allocate(coarse_state(nvar, n, n, n), coarse_temperature(n, n, n))
    allocate(fine_state(nvar, patch%fine_nx(), patch%fine_ny(), &
      patch%fine_nz()))
    allocate(fine_temperature(patch%fine_nx(), patch%fine_ny(), &
      patch%fine_nz()))
    allocate(x(n), y(n), z(n), xf(patch%fine_nx()), &
      yf(patch%fine_ny()), zf(patch%fine_nz()))
    allocate(mass_fractions(size(species)))
    call uniform_cell_centers_3d( &
      n, n, n, config%x_lower, config%x_upper, &
      config%y_lower, config%y_upper, config%z_lower, config%z_upper, &
      x, y, z, dx, dy, dz)
    call fine_patch_centers(patch, config, dx, dy, dz, xf, yf, zf)
    call initialize_reactive_problem_3d( &
      species, config, x, y, z, coarse_state, coarse_temperature, &
      density, mass_fractions, local_ok)
    if (.not. local_ok) return
    fine_config = config
    fine_config%nx = patch%fine_nx()
    fine_config%ny = patch%fine_ny()
    fine_config%nz = patch%fine_nz()
    call initialize_reactive_problem_3d( &
      species, fine_config, xf, yf, zf, fine_state, fine_temperature, &
      fine_density, mass_fractions, local_ok)
    if (.not. local_ok) return
    ok_out = density > 0.0_dp .and. fine_density > 0.0_dp .and. &
      all(ieee_is_finite(coarse_state)) .and. &
      all(ieee_is_finite(coarse_temperature)) .and. &
      all(ieee_is_finite(fine_state)) .and. &
      all(ieee_is_finite(fine_temperature))
  end subroutine initialize_public_api_hierarchy

  subroutine run_entropy_wave( &
      species, patch, reconstruction, conservation_error, synchronization_error, &
      maximum_reflux, closure_error, density_error, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    character(len=*), intent(in) :: reconstruction
    real(dp), intent(out) :: conservation_error, synchronization_error
    real(dp), intent(out) :: maximum_reflux, closure_error
    real(dp), intent(out) :: density_error
    logical, intent(out) :: ok

    type(reactive_3d_config) :: config
    real(dp), allocatable :: coarse_state(:, :, :, :)
    real(dp), allocatable :: fine_state(:, :, :, :)
    real(dp), allocatable :: saved_coarse_state(:, :, :, :)
    real(dp), allocatable :: saved_fine_state(:, :, :, :)
    real(dp), allocatable :: coarse_temperature(:, :, :)
    real(dp), allocatable :: fine_temperature(:, :, :)
    real(dp), allocatable :: saved_coarse_temperature(:, :, :)
    real(dp), allocatable :: saved_fine_temperature(:, :, :)
    real(dp), allocatable :: recovered_temperature(:, :, :)
    real(dp), allocatable :: restricted(:, :, :, :)
    real(dp), allocatable :: initial_integrals(:), final_integrals(:)
    real(dp), allocatable :: mass_fractions(:)
    real(dp), allocatable :: x(:), y(:), z(:), xf(:), yf(:), zf(:)
    real(dp) :: dx, dy, dz
    real(dp) :: coarse_density, fine_density, time, dt, reflux
    real(dp) :: coarse_closure, fine_closure, state_change
    real(dp) :: exact_density
    logical :: local_ok
    integer :: nvar, steps, i, j, k

    config = reactive_3d_config()
    config%nx = n
    config%ny = n
    config%nz = n
    config%thermo_model = "elementary"
    config%riemann_solver = "pelec"
    config%initial_velocity_x = 240.0_dp
    config%initial_velocity_y = -90.0_dp
    config%initial_velocity_z = 60.0_dp
    config%wave_number_x = 1
    config%wave_number_y = 2
    config%wave_number_z = 1
    nvar = reactive_nvar(size(species))
    allocate(coarse_state(nvar, n, n, n), coarse_temperature(n, n, n))
    allocate(fine_state(nvar, patch%fine_nx(), patch%fine_ny(), &
      patch%fine_nz()))
    allocate(fine_temperature(patch%fine_nx(), patch%fine_ny(), &
      patch%fine_nz()))
    allocate(x(n), y(n), z(n), xf(patch%fine_nx()), &
      yf(patch%fine_ny()), zf(patch%fine_nz()))
    allocate(mass_fractions(size(species)))
    allocate(initial_integrals(nvar), final_integrals(nvar))
    call uniform_cell_centers_3d( &
      n, n, n, config%x_lower, config%x_upper, &
      config%y_lower, config%y_upper, config%z_lower, config%z_upper, &
      x, y, z, dx, dy, dz)
    call fine_patch_centers(patch, config, dx, dy, dz, xf, yf, zf)
    call initialize_reactive_problem_3d( &
      species, config, x, y, z, coarse_state, coarse_temperature, &
      coarse_density, mass_fractions, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    call initialize_reactive_problem_3d( &
      species, config, xf, yf, zf, fine_state, fine_temperature, &
      fine_density, mass_fractions, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    call average_down_3d(coarse_state, fine_state, patch, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    allocate(recovered_temperature(n, n, n))
    call recover_reactive_temperatures_3d( &
      species, coarse_state, coarse_temperature, n, n, n, &
      recovered_temperature, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    coarse_temperature = recovered_temperature
    call composite_integrals_amr_3d( &
      coarse_state, fine_state, patch, dx, dy, dz, initial_integrals, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    saved_coarse_state = coarse_state
    maximum_reflux = 0.0_dp
    time = 0.0_dp
    steps = 0
    do while (time < final_time - 1.0e-15_dp)
      call compute_amr_reactive_cfl_timestep_3d( &
        species, patch, coarse_state, coarse_temperature, &
        fine_state, fine_temperature, dx, dy, dz, 0.30_dp, dt, local_ok)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
      dt = min(dt, final_time - time)
      call advance_amr_reactive_hydro_3d( &
        species, patch, coarse_state, coarse_temperature, &
        fine_state, fine_temperature, dx, dy, dz, dt, &
        config%riemann_solver, reflux, local_ok, reconstruction, "mc")
      if (.not. local_ok) then
        ok = .false.
        return
      end if
      maximum_reflux = max(maximum_reflux, reflux)
      time = time + dt
      steps = steps + 1
      if (steps > 10000) then
        ok = .false.
        return
      end if
    end do
    state_change = maxval(abs(coarse_state - saved_coarse_state))
    call composite_integrals_amr_3d( &
      coarse_state, fine_state, patch, dx, dy, dz, final_integrals, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    conservation_error = maxval(abs(final_integrals - initial_integrals) / &
      max(1.0_dp, abs(initial_integrals)))
    allocate(restricted(nvar, 4, 4, 4))
    call restrict_average_3d(fine_state, patch, restricted, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    synchronization_error = maxval(abs(restricted - &
      coarse_state(:, 3:6, 3:6, 3:6))) / &
      max(1.0_dp, maxval(abs(restricted)))
    call extrema_closure( &
      species, coarse_state, coarse_temperature, coarse_closure, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    call extrema_closure( &
      species, fine_state, fine_temperature, fine_closure, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    closure_error = max(coarse_closure, fine_closure)
    density_error = 0.0_dp
    do k = 1, size(zf)
      do j = 1, size(yf)
        do i = 1, size(xf)
          exact_density = reactive_entropy_wave_density_3d( &
            xf(i), yf(j), zf(k), time, config, coarse_density)
          density_error = density_error + &
            abs(fine_state(1, i, j, k) - exact_density)
        end do
      end do
    end do
    density_error = density_error / &
      real(size(xf) * size(yf) * size(zf), dp)

    saved_coarse_state = coarse_state
    saved_coarse_temperature = coarse_temperature
    saved_fine_state = fine_state
    saved_fine_temperature = fine_temperature
    call advance_amr_reactive_hydro_3d( &
      species, patch, coarse_state, coarse_temperature, &
      fine_state, fine_temperature, dx, dy, dz, 1.0e-9_dp, &
      "invalid", reflux, local_ok, reconstruction, "mc")
    ok = .not. local_ok .and. reflux == 0.0_dp .and. &
      maxval(abs(coarse_state - saved_coarse_state)) == 0.0_dp .and. &
      maxval(abs(coarse_temperature - saved_coarse_temperature)) == 0.0_dp &
      .and. maxval(abs(fine_state - saved_fine_state)) == 0.0_dp .and. &
      maxval(abs(fine_temperature - saved_fine_temperature)) == 0.0_dp &
      .and. state_change > 0.0_dp .and. coarse_density > 0.0_dp .and. &
      fine_density > 0.0_dp
    if (.not. ok) return
    call advance_amr_reactive_hydro_3d( &
      species, patch, coarse_state, coarse_temperature, &
      fine_state, fine_temperature, dx, dy, dz, 1.0e-9_dp, &
      "pelec", reflux, local_ok, "invalid", "mc")
    ok = .not. local_ok .and. reflux == 0.0_dp .and. &
      maxval(abs(coarse_state - saved_coarse_state)) == 0.0_dp .and. &
      maxval(abs(coarse_temperature - saved_coarse_temperature)) == 0.0_dp &
      .and. maxval(abs(fine_state - saved_fine_state)) == 0.0_dp .and. &
      maxval(abs(fine_temperature - saved_fine_temperature)) == 0.0_dp
    if (.not. ok) return
    call advance_amr_reactive_hydro_3d( &
      species, patch, coarse_state, coarse_temperature, &
      fine_state, fine_temperature, dx, dy, dz, 1.0e-9_dp, &
      "pelec", reflux, local_ok, "characteristic_plm", "invalid")
    ok = .not. local_ok .and. reflux == 0.0_dp .and. &
      maxval(abs(coarse_state - saved_coarse_state)) == 0.0_dp .and. &
      maxval(abs(coarse_temperature - saved_coarse_temperature)) == 0.0_dp &
      .and. maxval(abs(fine_state - saved_fine_state)) == 0.0_dp .and. &
      maxval(abs(fine_temperature - saved_fine_temperature)) == 0.0_dp
  end subroutine run_entropy_wave

  subroutine run_uniform_invariance( &
      species, patch, reconstruction, limiter, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    character(len=*), intent(in) :: reconstruction
    character(len=*), intent(in) :: limiter
    logical, intent(out) :: ok

    type(reactive_3d_config) :: config
    real(dp), allocatable :: coarse_state(:, :, :, :), fine_state(:, :, :, :)
    real(dp), allocatable :: initial_coarse(:, :, :, :), initial_fine(:, :, :, :)
    real(dp), allocatable :: coarse_temperature(:, :, :)
    real(dp), allocatable :: fine_temperature(:, :, :)
    real(dp), allocatable :: x(:), y(:), z(:), xf(:), yf(:), zf(:)
    real(dp), allocatable :: mass_fractions(:)
    real(dp) :: dx, dy, dz, density, fine_density, reflux, dt
    logical :: local_ok
    integer :: nvar

    config = reactive_3d_config()
    config%nx = n
    config%ny = n
    config%nz = n
    config%thermo_model = "elementary"
    config%density_wave_amplitude = 0.0_dp
    nvar = reactive_nvar(size(species))
    allocate(coarse_state(nvar, n, n, n), coarse_temperature(n, n, n))
    allocate(fine_state(nvar, patch%fine_nx(), patch%fine_ny(), &
      patch%fine_nz()))
    allocate(fine_temperature(patch%fine_nx(), patch%fine_ny(), &
      patch%fine_nz()))
    allocate(x(n), y(n), z(n), xf(patch%fine_nx()), &
      yf(patch%fine_ny()), zf(patch%fine_nz()))
    allocate(mass_fractions(size(species)))
    call uniform_cell_centers_3d( &
      n, n, n, config%x_lower, config%x_upper, &
      config%y_lower, config%y_upper, config%z_lower, config%z_upper, &
      x, y, z, dx, dy, dz)
    call fine_patch_centers(patch, config, dx, dy, dz, xf, yf, zf)
    call initialize_reactive_problem_3d( &
      species, config, x, y, z, coarse_state, coarse_temperature, &
      density, mass_fractions, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    call initialize_reactive_problem_3d( &
      species, config, xf, yf, zf, fine_state, fine_temperature, &
      fine_density, mass_fractions, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    initial_coarse = coarse_state
    initial_fine = fine_state
    call compute_amr_reactive_cfl_timestep_3d( &
      species, patch, coarse_state, coarse_temperature, &
      fine_state, fine_temperature, dx, dy, dz, 0.30_dp, dt, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    call advance_amr_reactive_hydro_3d( &
      species, patch, coarse_state, coarse_temperature, &
      fine_state, fine_temperature, dx, dy, dz, dt, "rusanov", &
      reflux, local_ok, reconstruction, limiter)
    ok = local_ok .and. &
      maxval(abs(coarse_state - initial_coarse)) <= &
        5.0e-14_dp * max(1.0_dp, maxval(abs(initial_coarse))) .and. &
      maxval(abs(fine_state - initial_fine)) <= &
        5.0e-14_dp * max(1.0_dp, maxval(abs(initial_fine))) .and. &
      reflux <= 5.0e-14_dp * max(1.0_dp, maxval(abs(initial_coarse))) .and. &
      density > 0.0_dp .and. fine_density > 0.0_dp
  end subroutine run_uniform_invariance

  subroutine fine_patch_centers(patch, config, dx, dy, dz, x, y, z)
    type(amr_patch_3d), intent(in) :: patch
    type(reactive_3d_config), intent(in) :: config
    real(dp), intent(in) :: dx, dy, dz
    real(dp), intent(out) :: x(:), y(:), z(:)
    integer :: i

    do i = 1, size(x)
      x(i) = config%x_lower + real(patch%coarse_i_lower - 1, dp) * dx + &
        (real(i, dp) - 0.5_dp) * dx / real(patch%refinement_ratio, dp)
    end do
    do i = 1, size(y)
      y(i) = config%y_lower + real(patch%coarse_j_lower - 1, dp) * dy + &
        (real(i, dp) - 0.5_dp) * dy / real(patch%refinement_ratio, dp)
    end do
    do i = 1, size(z)
      z(i) = config%z_lower + real(patch%coarse_k_lower - 1, dp) * dz + &
        (real(i, dp) - 0.5_dp) * dz / real(patch%refinement_ratio, dp)
    end do
  end subroutine fine_patch_centers

  subroutine extrema_closure(species, state, temperature, closure, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    real(dp), intent(out) :: closure
    logical, intent(out) :: ok
    real(dp) :: minimum_density, maximum_density
    real(dp) :: minimum_pressure, maximum_pressure
    real(dp) :: minimum_temperature, maximum_temperature, maximum_speed

    call reactive_extrema_3d( &
      species, state, temperature, size(state, 2), size(state, 3), &
      size(state, 4), minimum_density, maximum_density, minimum_pressure, &
      maximum_pressure, minimum_temperature, maximum_temperature, &
      maximum_speed, closure, ok)
    ok = ok .and. minimum_density > 0.0_dp .and. &
      minimum_pressure > 0.0_dp .and. minimum_temperature > 0.0_dp .and. &
      maximum_density > 0.0_dp .and. maximum_pressure > 0.0_dp .and. &
      maximum_temperature > 0.0_dp .and. maximum_speed >= 0.0_dp
  end subroutine extrema_closure

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) then
      write(*, '(a)') "FAIL: " // trim(message)
      error stop 1
    end if
  end subroutine require

end program test_amr_reactive_3d
