program test_reactive_euler_3d_dimensional_reduction
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use state_indices_mod, only: imx, imy, imz
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: &
    mass_fractions_from_mole_fractions, mixture_density
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_primitive_to_conserved, reactive_conserved_to_primitive, &
    reactive_riemann_flux_x
  use reactive_3d_mod, only: &
    advance_reactive_euler_ssprk2_3d, &
    advance_reactive_euler_ssprk2_with_fluxes_3d, &
    advance_reactive_euler_ssprk2_plm_3d, &
    advance_reactive_euler_ssprk2_plm_with_fluxes_3d
  use simulation_config_reactive_3d_mod, only: &
    reactive_3d_config, validate_reactive_3d_configuration, &
    reactive_3d_mole_fractions
  implicit none

  integer, parameter :: n = 6
  real(dp), parameter :: dx = 0.01_dp / real(n, dp)
  real(dp), parameter :: dt = 2.0e-7_dp
  type(nasa7_species), allocatable :: species(:)
  real(dp), allocatable :: base_state(:, :), base_temperature(:)
  real(dp), allocatable :: reference_state(:, :), reference_temperature(:)
  character(len=8), parameter :: solvers(3) = [ &
    character(len=8) :: "rusanov", "hllc", "pelec"]
  real(dp) :: maximum_state_error, maximum_temperature_error
  logical :: ok
  integer :: solver_index, direction

  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok) error stop "Failed to load elementary thermodynamics"
  call check_configuration_contract(ok)
  if (.not. ok) error stop "Reactive 3D configuration contract failed"
  call build_line_state(species, base_state, base_temperature, ok)
  if (.not. ok) error stop "Failed to build dimensional-reduction state"

  do solver_index = 1, size(solvers)
    reference_state = base_state
    reference_temperature = base_temperature
    call advance_reference_ssprk2_x( &
      species, reference_state, reference_temperature, dx, dt, &
      trim(solvers(solver_index)), ok)
    if (.not. ok) error stop "Independent 1D SSPRK2 reference failed"
    do direction = 1, 3
      call compare_direction( &
        species, base_state, base_temperature, reference_state, &
        reference_temperature, direction, trim(solvers(solver_index)), &
        maximum_state_error, maximum_temperature_error, ok)
      if (.not. ok) error stop "Reactive 3D dimensional reduction failed"
      write(*, '(a,1x,a,a,i0,2(a,es24.16))') &
        "solver", trim(solvers(solver_index)), ", direction=", direction, &
        ", state error=", maximum_state_error, &
        ", temperature error=", maximum_temperature_error
      if (maximum_state_error > 5.0e-13_dp .or. &
          maximum_temperature_error > 5.0e-13_dp) then
        error stop "Reactive 3D dimensional-reduction tolerance exceeded"
      end if
    end do
  end do

  call check_transactional_rejection( &
    species, base_state, base_temperature, ok)
  if (.not. ok) error stop "Reactive 3D rejection was not transactional"
  call check_flux_returning_equivalence( &
    species, base_state, base_temperature, ok)
  if (.not. ok) error stop "Reactive 3D flux-returning advance mismatch"
  call check_plm_directional_symmetry( &
    species, base_state, base_temperature, ok)
  if (.not. ok) error stop "Reactive 3D PLM directional symmetry failed"
  call check_plm_transactional_rejection( &
    species, base_state, base_temperature, ok)
  if (.not. ok) error stop "Reactive 3D PLM rejection was not transactional"
  write(*, '(a)') "test_reactive_euler_3d_dimensional_reduction: PASS"

contains

  subroutine check_configuration_contract(ok)
    logical, intent(out) :: ok

    type(reactive_3d_config) :: config
    real(dp) :: full_fractions(10), elementary_fractions(7)
    character(len=256) :: message
    logical :: local_ok

    config = reactive_3d_config()
    call validate_reactive_3d_configuration(config, local_ok, message)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    config%reconstruction = "characteristic_plm"
    config%limiter = "mc"
    call validate_reactive_3d_configuration(config, local_ok, message)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    config%nx = 2
    call validate_reactive_3d_configuration(config, local_ok, message)
    if (local_ok) then
      ok = .false.
      return
    end if
    config%nx = 12
    config%reconstruction = "invalid"
    call validate_reactive_3d_configuration(config, local_ok, message)
    if (local_ok) then
      ok = .false.
      return
    end if
    config%reconstruction = "characteristic_plm"
    config%limiter = "invalid"
    call validate_reactive_3d_configuration(config, local_ok, message)
    if (local_ok) then
      ok = .false.
      return
    end if
    config%reconstruction = "pcm"
    config%limiter = "mc"
    call reactive_3d_mole_fractions(config, 10, full_fractions, local_ok)
    if (.not. local_ok .or. abs(sum(full_fractions) - 1.0_dp) > 5.0e-14_dp) then
      ok = .false.
      return
    end if
    config%thermo_model = "elementary"
    call validate_reactive_3d_configuration(config, local_ok, message)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    call reactive_3d_mole_fractions( &
      config, 7, elementary_fractions, local_ok)
    if (.not. local_ok .or. &
        abs(sum(elementary_fractions) - 1.0_dp) > 5.0e-14_dp) then
      ok = .false.
      return
    end if
    config%boundary_condition = "outflow"
    call validate_reactive_3d_configuration(config, local_ok, message)
    if (local_ok) then
      ok = .false.
      return
    end if
    config%boundary_condition = "periodic"
    config%riemann_solver = "invalid"
    call validate_reactive_3d_configuration(config, local_ok, message)
    if (local_ok) then
      ok = .false.
      return
    end if
    config%riemann_solver = "rusanov"
    config%x_n2 = config%x_n2 + 0.01_dp
    call validate_reactive_3d_configuration(config, local_ok, message)
    if (local_ok) then
      ok = .false.
      return
    end if
    config = reactive_3d_config()
    config%problem = "reactive_hotspot"
    config%chemistry_enabled = .true.
    config%hotspot_width = 0.0_dp
    call validate_reactive_3d_configuration(config, local_ok, message)
    if (local_ok) then
      ok = .false.
      return
    end if
    config%hotspot_width = 0.0012_dp
    config%chemistry_relative_tolerance = 0.0_dp
    call validate_reactive_3d_configuration(config, local_ok, message)
    if (local_ok) then
      ok = .false.
      return
    end if
    config%chemistry_relative_tolerance = 2.0e-7_dp
    config%transport_enabled = .true.
    config%viscosity_enabled = .false.
    config%thermal_conduction_enabled = .false.
    config%species_diffusion_enabled = .false.
    config%barodiffusion_enabled = .false.
    call validate_reactive_3d_configuration(config, local_ok, message)
    if (local_ok) then
      ok = .false.
      return
    end if
    config%thermal_conduction_enabled = .true.
    config%transport_cfl = 0.5_dp
    call validate_reactive_3d_configuration(config, local_ok, message)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    config%transport_cfl = 0.500001_dp
    call validate_reactive_3d_configuration(config, local_ok, message)
    if (local_ok) then
      ok = .false.
      return
    end if
    config%transport_cfl = 0.35_dp
    call validate_reactive_3d_configuration(config, ok, message)
  end subroutine check_configuration_contract

  subroutine build_line_state(species, state, temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), allocatable, intent(out) :: state(:, :), temperature(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: mole_fractions(:), mass_fractions(:), primitive(:)
    real(dp) :: base_density, density, sound_speed, phase
    logical :: cell_ok
    integer :: i, species_index, nvar

    nvar = reactive_nvar(size(species))
    allocate(state(nvar, n), temperature(n))
    allocate(mole_fractions(size(species)), mass_fractions(size(species)))
    allocate(primitive(reactive_nprim(size(species))))
    mole_fractions = [ &
      0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
      1.0e-5_dp, 0.0_dp, 0.55643_dp]
    call mass_fractions_from_mole_fractions( &
      species, mole_fractions, mass_fractions, cell_ok)
    if (.not. cell_ok) then
      ok = .false.
      return
    end if
    base_density = mixture_density( &
      species, mass_fractions, 101325.0_dp, 1000.0_dp, cell_ok)
    if (.not. cell_ok) then
      ok = .false.
      return
    end if

    do i = 1, n
      phase = 2.0_dp * acos(-1.0_dp) * &
        (real(i, dp) - 0.5_dp) / real(n, dp)
      density = base_density * (1.0_dp + 0.05_dp * sin(phase))
      primitive(1:5) = [density, 120.0_dp, -30.0_dp, 20.0_dp, 101325.0_dp]
      do species_index = 1, size(species)
        primitive(reactive_mass_fraction_component(species_index)) = &
          mass_fractions(species_index)
      end do
      call reactive_primitive_to_conserved( &
        species, primitive, state(:, i), temperature(i), sound_speed, cell_ok)
      if (.not. cell_ok) then
        ok = .false.
        return
      end if
    end do
    ok = .true.
  end subroutine build_line_state

  subroutine advance_reference_ssprk2_x( &
      species, state, temperature, spacing, timestep, solver, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:, :), temperature(:)
    real(dp), intent(in) :: spacing, timestep
    character(len=*), intent(in) :: solver
    logical, intent(out) :: ok

    real(dp), allocatable :: old_state(:, :), stage_state(:, :), rhs(:, :)
    real(dp), allocatable :: old_temperature(:), stage_temperature(:)
    real(dp), allocatable :: updated_temperature(:)
    logical :: local_ok
    integer :: nvar

    nvar = size(state, 1)
    allocate(old_state(nvar, n), stage_state(nvar, n), rhs(nvar, n))
    allocate(old_temperature(n), stage_temperature(n), updated_temperature(n))
    old_state = state
    call recover_line( &
      species, old_state, temperature, old_temperature, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    call reference_rhs_x( &
      species, old_state, old_temperature, spacing, solver, rhs, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    stage_state = old_state + timestep * rhs
    call recover_line( &
      species, stage_state, old_temperature, stage_temperature, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    call reference_rhs_x( &
      species, stage_state, stage_temperature, spacing, solver, rhs, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    state = 0.5_dp * old_state + 0.5_dp * (stage_state + timestep * rhs)
    call recover_line( &
      species, state, stage_temperature, updated_temperature, ok)
    if (ok) temperature = updated_temperature
  end subroutine advance_reference_ssprk2_x

  subroutine reference_rhs_x( &
      species, state, temperature, spacing, solver, rhs, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :), temperature(:), spacing
    character(len=*), intent(in) :: solver
    real(dp), intent(out) :: rhs(:, :)
    logical, intent(out) :: ok

    real(dp), allocatable :: flux(:, :)
    logical :: face_ok
    integer :: i, next_i, previous_i

    rhs = 0.0_dp
    allocate(flux(size(state, 1), n))
    do i = 1, n
      next_i = modulo(i, n) + 1
      call reactive_riemann_flux_x( &
        species, state(:, i), state(:, next_i), temperature(i), &
        temperature(next_i), solver, flux(:, i), face_ok)
      if (.not. face_ok) then
        ok = .false.
        return
      end if
    end do
    do i = 1, n
      previous_i = modulo(i - 2, n) + 1
      rhs(:, i) = -(flux(:, i) - flux(:, previous_i)) / spacing
    end do
    ok = .true.
  end subroutine reference_rhs_x

  subroutine recover_line(species, state, guess, temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :), guess(:)
    real(dp), intent(out) :: temperature(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:)
    real(dp) :: sound_speed
    logical :: cell_ok
    integer :: i

    allocate(primitive(reactive_nprim(size(species))))
    do i = 1, n
      call reactive_conserved_to_primitive( &
        species, state(:, i), guess(i), primitive, temperature(i), &
        sound_speed, cell_ok)
      if (.not. cell_ok) then
        ok = .false.
        return
      end if
    end do
    ok = .true.
  end subroutine recover_line

  subroutine compare_direction( &
      species, base_state, base_temperature, reference_state, &
      reference_temperature, direction, solver, state_error, &
      temperature_error, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: base_state(:, :), base_temperature(:)
    real(dp), intent(in) :: reference_state(:, :), reference_temperature(:)
    integer, intent(in) :: direction
    character(len=*), intent(in) :: solver
    real(dp), intent(out) :: state_error, temperature_error
    logical, intent(out) :: ok

    real(dp), allocatable :: grid_state(:, :, :, :), grid_temperature(:, :, :)
    real(dp), allocatable :: expected(:), rotated(:)
    real(dp) :: scale, local_state_error, local_temperature_error
    integer :: nx, ny, nz, i, j, k, line_index, nvar

    nvar = size(base_state, 1)
    select case (direction)
    case (1)
      nx = n
      ny = 2
      nz = 2
    case (2)
      nx = 2
      ny = n
      nz = 2
    case (3)
      nx = 2
      ny = 2
      nz = n
    case default
      ok = .false.
      return
    end select
    allocate(grid_state(nvar, nx, ny, nz), grid_temperature(nx, ny, nz))
    allocate(expected(nvar), rotated(nvar))
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          line_index = merge(i, merge(j, k, direction == 2), direction == 1)
          call rotate_from_x(base_state(:, line_index), direction, rotated)
          grid_state(:, i, j, k) = rotated
          grid_temperature(i, j, k) = base_temperature(line_index)
        end do
      end do
    end do

    call advance_reactive_euler_ssprk2_3d( &
      species, grid_state, grid_temperature, nx, ny, nz, &
      merge(dx, 1.3_dp * dx, direction == 1), &
      merge(dx, 1.1_dp * dx, direction == 2), &
      merge(dx, 0.9_dp * dx, direction == 3), dt, solver, ok)
    if (.not. ok) return

    state_error = 0.0_dp
    temperature_error = 0.0_dp
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          line_index = merge(i, merge(j, k, direction == 2), direction == 1)
          call rotate_from_x(reference_state(:, line_index), direction, expected)
          scale = max(1.0_dp, maxval(abs(expected)))
          local_state_error = &
            maxval(abs(grid_state(:, i, j, k) - expected)) / scale
          local_temperature_error = &
            abs(grid_temperature(i, j, k) - &
              reference_temperature(line_index)) / &
            max(1.0_dp, abs(reference_temperature(line_index)))
          state_error = max(state_error, local_state_error)
          temperature_error = max( &
            temperature_error, local_temperature_error)
        end do
      end do
    end do
  end subroutine compare_direction

  pure subroutine rotate_from_x(input, direction, output)
    real(dp), intent(in) :: input(:)
    integer, intent(in) :: direction
    real(dp), intent(out) :: output(:)

    output = input
    select case (direction)
    case (2)
      output(imx) = input(imy)
      output(imy) = input(imx)
    case (3)
      output(imx) = input(imz)
      output(imz) = input(imx)
    end select
  end subroutine rotate_from_x

  subroutine check_transactional_rejection( &
      species, base_state, base_temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: base_state(:, :), base_temperature(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: state(:, :, :, :), temperature(:, :, :)
    real(dp), allocatable :: saved_state(:, :, :, :), saved_temperature(:, :, :)
    integer :: i, j, k, nvar

    nvar = size(base_state, 1)
    allocate(state(nvar, n, 2, 2), temperature(n, 2, 2))
    do k = 1, 2
      do j = 1, 2
        do i = 1, n
          state(:, i, j, k) = base_state(:, i)
          temperature(i, j, k) = base_temperature(i)
        end do
      end do
    end do
    saved_state = state
    saved_temperature = temperature
    call advance_reactive_euler_ssprk2_3d( &
      species, state, temperature, n, 2, 2, dx, dx, dx, dt, &
      "invalid", ok)
    ok = .not. ok .and. maxval(abs(state - saved_state)) == 0.0_dp .and. &
      maxval(abs(temperature - saved_temperature)) == 0.0_dp
  end subroutine check_transactional_rejection

  subroutine check_flux_returning_equivalence( &
      species, base_state, base_temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: base_state(:, :), base_temperature(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: reference_state(:, :, :, :)
    real(dp), allocatable :: flux_state(:, :, :, :)
    real(dp), allocatable :: reference_temperature(:, :, :)
    real(dp), allocatable :: flux_temperature(:, :, :)
    real(dp), allocatable :: flux_x(:, :, :, :)
    real(dp), allocatable :: flux_y(:, :, :, :)
    real(dp), allocatable :: flux_z(:, :, :, :)
    real(dp), allocatable :: saved_state(:, :, :, :)
    real(dp), allocatable :: saved_temperature(:, :, :)
    logical :: local_ok
    integer :: i, j, k, nvar

    nvar = size(base_state, 1)
    allocate(reference_state(nvar, n, 2, 2))
    allocate(reference_temperature(n, 2, 2))
    do k = 1, 2
      do j = 1, 2
        do i = 1, n
          reference_state(:, i, j, k) = base_state(:, i)
          reference_temperature(i, j, k) = base_temperature(i)
        end do
      end do
    end do
    flux_state = reference_state
    flux_temperature = reference_temperature
    allocate(flux_x, mold=reference_state)
    allocate(flux_y, mold=reference_state)
    allocate(flux_z, mold=reference_state)
    call advance_reactive_euler_ssprk2_3d( &
      species, reference_state, reference_temperature, n, 2, 2, &
      dx, 1.1_dp * dx, 0.9_dp * dx, dt, "pelec", local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    call advance_reactive_euler_ssprk2_with_fluxes_3d( &
      species, flux_state, flux_temperature, n, 2, 2, &
      dx, 1.1_dp * dx, 0.9_dp * dx, dt, "pelec", &
      flux_x, flux_y, flux_z, local_ok)
    ok = local_ok .and. &
      maxval(abs(flux_state - reference_state)) == 0.0_dp .and. &
      maxval(abs(flux_temperature - reference_temperature)) == 0.0_dp .and. &
      all(ieee_is_finite(flux_x)) .and. all(ieee_is_finite(flux_y)) .and. &
      all(ieee_is_finite(flux_z))
    if (.not. ok) return

    saved_state = flux_state
    saved_temperature = flux_temperature
    call advance_reactive_euler_ssprk2_with_fluxes_3d( &
      species, flux_state, flux_temperature, n, 2, 2, &
      dx, dx, dx, dt, "invalid", flux_x, flux_y, flux_z, local_ok)
    ok = .not. local_ok .and. &
      maxval(abs(flux_state - saved_state)) == 0.0_dp .and. &
      maxval(abs(flux_temperature - saved_temperature)) == 0.0_dp .and. &
      maxval(abs(flux_x)) == 0.0_dp .and. &
      maxval(abs(flux_y)) == 0.0_dp .and. &
      maxval(abs(flux_z)) == 0.0_dp
  end subroutine check_flux_returning_equivalence

  subroutine check_plm_directional_symmetry( &
      species, base_state, base_temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: base_state(:, :), base_temperature(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: reference_state(:, :), candidate_state(:, :)
    real(dp), allocatable :: reference_temperature(:)
    real(dp), allocatable :: candidate_temperature(:), expected(:)
    real(dp) :: state_error, temperature_error, scale
    logical :: local_ok
    integer :: direction, i

    call run_plm_direction( &
      species, base_state, base_temperature, 1, reference_state, &
      reference_temperature, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    allocate(expected(size(base_state, 1)))
    do direction = 2, 3
      call run_plm_direction( &
        species, base_state, base_temperature, direction, candidate_state, &
        candidate_temperature, local_ok)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
      state_error = 0.0_dp
      temperature_error = 0.0_dp
      do i = 1, n
        call rotate_from_x(reference_state(:, i), direction, expected)
        scale = max(1.0_dp, maxval(abs(expected)))
        state_error = max( &
          state_error, maxval(abs(candidate_state(:, i) - expected)) / scale)
        temperature_error = max( &
          temperature_error, &
          abs(candidate_temperature(i) - reference_temperature(i)) / &
            max(1.0_dp, abs(reference_temperature(i))))
      end do
      write(*, '(a,i0,2(a,es24.16))') &
        "PLM direction=", direction, ", state error=", state_error, &
        ", temperature error=", temperature_error
      if (state_error > 5.0e-13_dp .or. &
          temperature_error > 5.0e-13_dp) then
        ok = .false.
        return
      end if
      deallocate(candidate_state, candidate_temperature)
    end do
    ok = .true.
  end subroutine check_plm_directional_symmetry

  subroutine run_plm_direction( &
      species, base_state, base_temperature, direction, line_state, &
      line_temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: base_state(:, :), base_temperature(:)
    integer, intent(in) :: direction
    real(dp), allocatable, intent(out) :: line_state(:, :)
    real(dp), allocatable, intent(out) :: line_temperature(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: grid_state(:, :, :, :), grid_temperature(:, :, :)
    real(dp), allocatable :: rotated(:)
    integer :: nx, ny, nz, i, j, k, line_index, nvar

    select case (direction)
    case (1)
      nx = n
      ny = 3
      nz = 3
    case (2)
      nx = 3
      ny = n
      nz = 3
    case (3)
      nx = 3
      ny = 3
      nz = n
    case default
      ok = .false.
      return
    end select
    nvar = size(base_state, 1)
    allocate(grid_state(nvar, nx, ny, nz), grid_temperature(nx, ny, nz))
    allocate(rotated(nvar), line_state(nvar, n), line_temperature(n))
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          line_index = merge(i, merge(j, k, direction == 2), direction == 1)
          call rotate_from_x(base_state(:, line_index), direction, rotated)
          grid_state(:, i, j, k) = rotated
          grid_temperature(i, j, k) = base_temperature(line_index)
        end do
      end do
    end do
    call advance_reactive_euler_ssprk2_plm_3d( &
      species, grid_state, grid_temperature, nx, ny, nz, &
      merge(dx, 1.3_dp * dx, direction == 1), &
      merge(dx, 1.1_dp * dx, direction == 2), &
      merge(dx, 0.9_dp * dx, direction == 3), dt, "mc", "pelec", ok)
    if (.not. ok) return
    do line_index = 1, n
      select case (direction)
      case (1)
        line_state(:, line_index) = grid_state(:, line_index, 2, 2)
        line_temperature(line_index) = grid_temperature(line_index, 2, 2)
      case (2)
        line_state(:, line_index) = grid_state(:, 2, line_index, 2)
        line_temperature(line_index) = grid_temperature(2, line_index, 2)
      case (3)
        line_state(:, line_index) = grid_state(:, 2, 2, line_index)
        line_temperature(line_index) = grid_temperature(2, 2, line_index)
      end select
    end do
  end subroutine run_plm_direction

  subroutine check_plm_transactional_rejection( &
      species, base_state, base_temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: base_state(:, :), base_temperature(:)
    logical, intent(out) :: ok

    real(dp), allocatable :: state(:, :, :, :), temperature(:, :, :)
    real(dp), allocatable :: saved_state(:, :, :, :), saved_temperature(:, :, :)
    real(dp), allocatable :: flux_x(:, :, :, :), flux_y(:, :, :, :)
    real(dp), allocatable :: flux_z(:, :, :, :)
    logical :: local_ok
    integer :: i, j, k, nvar

    nvar = size(base_state, 1)
    allocate(state(nvar, n, 3, 3), temperature(n, 3, 3))
    allocate(flux_x(nvar, n, 3, 3), flux_y(nvar, n, 3, 3))
    allocate(flux_z(nvar, n, 3, 3))
    do k = 1, 3
      do j = 1, 3
        do i = 1, n
          state(:, i, j, k) = base_state(:, i)
          temperature(i, j, k) = base_temperature(i)
        end do
      end do
    end do
    saved_state = state
    saved_temperature = temperature
    call advance_reactive_euler_ssprk2_plm_with_fluxes_3d( &
      species, state, temperature, n, 3, 3, dx, dx, dx, dt, &
      "invalid", "pelec", flux_x, flux_y, flux_z, local_ok)
    ok = .not. local_ok .and. &
      maxval(abs(state - saved_state)) == 0.0_dp .and. &
      maxval(abs(temperature - saved_temperature)) == 0.0_dp .and. &
      maxval(abs(flux_x)) == 0.0_dp .and. &
      maxval(abs(flux_y)) == 0.0_dp .and. &
      maxval(abs(flux_z)) == 0.0_dp
    if (.not. ok) return
    call advance_reactive_euler_ssprk2_plm_with_fluxes_3d( &
      species, state, temperature, n, 3, 3, dx, dx, dx, dt, &
      "mc", "invalid", flux_x, flux_y, flux_z, local_ok)
    ok = .not. local_ok .and. &
      maxval(abs(state - saved_state)) == 0.0_dp .and. &
      maxval(abs(temperature - saved_temperature)) == 0.0_dp .and. &
      maxval(abs(flux_x)) == 0.0_dp .and. &
      maxval(abs(flux_y)) == 0.0_dp .and. &
      maxval(abs(flux_z)) == 0.0_dp
  end subroutine check_plm_transactional_rejection

end program test_reactive_euler_3d_dimensional_reduction
