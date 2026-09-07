program test_reactive_eb_hydro_3d
  use precision_mod, only: dp
  use state_indices_mod, only: irho
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_species_component, reactive_primitive_to_conserved
  use eb_geometry_3d_mod, only: &
    eb_geometry_3d, eb_covered_cell_3d, build_axis_plane_eb_geometry_3d
  use eb_reactive_wall_flux_3d_mod, only: &
    reactive_eb_flux_divergence_3d
  use eb_reactive_hydro_3d_mod, only: &
    reactive_eb_outflow_riemann_fluxes_3d, &
    advance_reactive_eb_euler_3d, &
    advance_reactive_eb_redistributed_euler_3d, &
    advance_reactive_eb_state_redistributed_euler_3d
  implicit none

  integer, parameter :: nx = 10
  integer, parameter :: ny = 8
  integer, parameter :: nz = 6
  real(dp), parameter :: pressure = 135000.0_dp
  type(nasa7_species), allocatable :: species(:)
  type(eb_geometry_3d) :: geometry
  real(dp), allocatable :: primitive(:), state_cell(:), perturbed_cell(:)
  real(dp), allocatable :: mass_fractions(:)
  real(dp), allocatable :: state(:, :, :, :), new_state(:, :, :, :)
  real(dp), allocatable :: temperature(:, :, :), new_temperature(:, :, :)
  real(dp), allocatable :: x_flux(:, :, :, :), y_flux(:, :, :, :)
  real(dp), allocatable :: z_flux(:, :, :, :)
  real(dp), allocatable :: before_integrals(:), after_integrals(:)
  real(dp) :: mole_fractions(7), base_temperature, perturbed_temperature
  real(dp) :: sound_speed, perturbed_sound_speed
  real(dp) :: maximum_uniform_state_error
  real(dp) :: maximum_uniform_temperature_error
  real(dp) :: conservation_error, species_closure_error
  real(dp) :: small_cell_full_step, redistributed_cut_density
  real(dp) :: state_redistributed_cut_density
  logical :: ok
  integer :: i, j, k, nvar, first_species, last_species

  maximum_uniform_state_error = 0.0_dp
  maximum_uniform_temperature_error = 0.0_dp

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "thermodynamic database load")
  nvar = reactive_nvar(size(species))
  first_species = reactive_species_component(1)
  last_species = reactive_species_component(size(species))
  allocate(primitive(reactive_nprim(size(species))))
  allocate(state_cell(nvar), perturbed_cell(nvar), mass_fractions(size(species)))
  allocate(state(nvar, nx, ny, nz), new_state(nvar, nx, ny, nz))
  allocate(temperature(nx, ny, nz), new_temperature(nx, ny, nz))
  allocate(x_flux(nvar, 0:nx, ny, nz))
  allocate(y_flux(nvar, nx, 0:ny, nz))
  allocate(z_flux(nvar, nx, ny, 0:nz))
  allocate(before_integrals(nvar), after_integrals(nvar))

  mole_fractions = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
    1.0e-5_dp, 0.0_dp, 0.55643_dp]
  call mass_fractions_from_mole_fractions( &
    species, mole_fractions, mass_fractions, ok)
  call require(ok, "composition conversion")

  call verify_uniform_tangent_update("x", "hllc", &
    [0.0_dp, 18.0_dp, -11.0_dp])
  call verify_uniform_tangent_update("y", "rusanov", &
    [18.0_dp, 0.0_dp, -11.0_dp])
  call verify_uniform_tangent_update("z", "pelec", &
    [18.0_dp, -11.0_dp, 0.0_dp])
  call verify_small_cell_hydro_relief()

  call build_axis_plane_eb_geometry_3d( &
    nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    "x", 1.1_dp, geometry, ok)
  call require(ok, "covered solver-validation geometry")
  call reactive_eb_outflow_riemann_fluxes_3d( &
    species, state, temperature, geometry, "unknown", &
    x_flux, y_flux, z_flux, ok)
  call require(.not. ok .and. maxval(abs(x_flux)) == 0.0_dp .and. &
    maxval(abs(y_flux)) == 0.0_dp .and. maxval(abs(z_flux)) == 0.0_dp, &
    "covered geometry invalid solver rejection")

  call build_axis_plane_eb_geometry_3d( &
    nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    "x", 0.37_dp, geometry, ok)
  call require(ok, "perturbation geometry")
  call build_cell_state(0.31_dp, [0.0_dp, 0.0_dp, 0.0_dp], &
    state_cell, base_temperature, sound_speed)
  call fill_state(state_cell, base_temperature)
  call build_cell_state(0.34_dp, [0.0_dp, 0.0_dp, 0.0_dp], &
    perturbed_cell, perturbed_temperature, perturbed_sound_speed)
  state(:, 7, 4, 3) = perturbed_cell
  temperature(7, 4, 3) = perturbed_temperature

  call reactive_eb_outflow_riemann_fluxes_3d( &
    species, state, temperature, geometry, "rusanov", &
    x_flux, y_flux, z_flux, ok)
  call require(ok, "perturbation face fluxes")
  call require(maxval(abs(x_flux(:, 3, :, :))) == 0.0_dp, &
    "closed embedded-normal face flux")
  call require(maxval(abs(x_flux(:, 4, :, :))) > 0.0_dp, &
    "open embedded-normal face flux")

  call fluid_integrals(state, geometry, before_integrals)
  call advance_reactive_eb_euler_3d( &
    species, state, temperature, geometry, "rusanov", 1.0e-8_dp, &
    new_state, new_temperature, ok)
  call require(ok, "perturbation Euler update")
  call fluid_integrals(new_state, geometry, after_integrals)
  conservation_error = maxval(abs(after_integrals - before_integrals) / &
    max(1.0_dp, abs(before_integrals)))
  call require(conservation_error <= 3.0e-13_dp, &
    "fluid-volume conservation")
  species_closure_error = maxval(abs( &
      sum(new_state(first_species:last_species, :, :, :), dim=1) - &
        new_state(irho, :, :, :)), &
      mask=geometry%cell_type /= eb_covered_cell_3d)
  call require(species_closure_error <= 3.0e-13_dp, &
    "species-density closure")
  call require(minval(new_temperature, &
    mask=geometry%cell_type /= eb_covered_cell_3d) > 0.0_dp, &
    "positive recovered temperature")
  call require(all(new_state(:, 1:3, :, :) == state(:, 1:3, :, :)) .and. &
    all(new_temperature(1:3, :, :) == temperature(1:3, :, :)), &
    "covered cells unchanged")

  call advance_reactive_eb_redistributed_euler_3d( &
    species, state, temperature, geometry, "rusanov", 1.0e-8_dp, &
    new_state, new_temperature, ok)
  call require(ok, "perturbation redistributed Euler update")
  call fluid_integrals(new_state, geometry, after_integrals)
  conservation_error = max(conservation_error, &
    maxval(abs(after_integrals - before_integrals) / &
      max(1.0_dp, abs(before_integrals))))
  call require(conservation_error <= 3.0e-13_dp, &
    "redistributed fluid-volume conservation")
  species_closure_error = max(species_closure_error, maxval(abs( &
      sum(new_state(first_species:last_species, :, :, :), dim=1) - &
        new_state(irho, :, :, :)), &
      mask=geometry%cell_type /= eb_covered_cell_3d))
  call require(species_closure_error <= 3.0e-13_dp, &
    "redistributed species-density closure")
  call require(minval(new_temperature, &
    mask=geometry%cell_type /= eb_covered_cell_3d) > 0.0_dp, &
    "positive redistributed temperature")
  call require(all(new_state(:, 1:3, :, :) == state(:, 1:3, :, :)) .and. &
    all(new_temperature(1:3, :, :) == temperature(1:3, :, :)), &
    "redistributed covered cells unchanged")

  call advance_reactive_eb_state_redistributed_euler_3d( &
    species, state, temperature, geometry, "rusanov", 1.0e-8_dp, &
    new_state, new_temperature, ok)
  call require(ok, "perturbation StateRedist Euler update")
  call fluid_integrals(new_state, geometry, after_integrals)
  conservation_error = max(conservation_error, &
    maxval(abs(after_integrals - before_integrals) / &
      max(1.0_dp, abs(before_integrals))))
  call require(conservation_error <= 3.0e-13_dp, &
    "StateRedist fluid-volume conservation")
  species_closure_error = max(species_closure_error, maxval(abs( &
      sum(new_state(first_species:last_species, :, :, :), dim=1) - &
        new_state(irho, :, :, :)), &
      mask=geometry%cell_type /= eb_covered_cell_3d))
  call require(species_closure_error <= 3.0e-13_dp, &
    "StateRedist species-density closure")
  call require(minval(new_temperature, &
    mask=geometry%cell_type /= eb_covered_cell_3d) > 0.0_dp, &
    "positive StateRedist temperature")
  call require(all(new_state(:, 1:3, :, :) == state(:, 1:3, :, :)) .and. &
    all(new_temperature(1:3, :, :) == temperature(1:3, :, :)), &
    "StateRedist covered cells unchanged")

  call advance_reactive_eb_euler_3d( &
    species, state, temperature, geometry, "unknown", 1.0e-8_dp, &
    new_state, new_temperature, ok)
  call require(.not. ok .and. all(new_state == state) .and. &
    all(new_temperature == temperature), "invalid solver rollback")
  call advance_reactive_eb_euler_3d( &
    species, state, temperature, geometry, "rusanov", -1.0e-8_dp, &
    new_state, new_temperature, ok)
  call require(.not. ok .and. all(new_state == state) .and. &
    all(new_temperature == temperature), "invalid timestep rollback")
  call advance_reactive_eb_redistributed_euler_3d( &
    species, state, temperature, geometry, "unknown", 1.0e-8_dp, &
    new_state, new_temperature, ok)
  call require(.not. ok .and. all(new_state == state) .and. &
    all(new_temperature == temperature), &
    "redistributed invalid solver rollback")
  call advance_reactive_eb_redistributed_euler_3d( &
    species, state, temperature, geometry, "rusanov", -1.0e-8_dp, &
    new_state, new_temperature, ok)
  call require(.not. ok .and. all(new_state == state) .and. &
    all(new_temperature == temperature), &
    "redistributed invalid timestep rollback")
  call advance_reactive_eb_state_redistributed_euler_3d( &
    species, state, temperature, geometry, "unknown", 1.0e-8_dp, &
    new_state, new_temperature, ok)
  call require(.not. ok .and. all(new_state == state) .and. &
    all(new_temperature == temperature), &
    "StateRedist invalid solver rollback")
  call advance_reactive_eb_state_redistributed_euler_3d( &
    species, state, temperature, geometry, "rusanov", -1.0e-8_dp, &
    new_state, new_temperature, ok)
  call require(.not. ok .and. all(new_state == state) .and. &
    all(new_temperature == temperature), &
    "StateRedist invalid timestep rollback")
  call advance_reactive_eb_state_redistributed_euler_3d( &
    species, state, temperature, geometry, "rusanov", 1.0e-8_dp, &
    new_state, new_temperature, ok, 0.0_dp)
  call require(.not. ok .and. all(new_state == state) .and. &
    all(new_temperature == temperature), &
    "StateRedist invalid target rollback")

  state(irho, 4, 2, 3) = -1.0_dp
  call advance_reactive_eb_euler_3d( &
    species, state, temperature, geometry, "rusanov", 1.0e-8_dp, &
    new_state, new_temperature, ok)
  call require(.not. ok .and. all(new_state == state) .and. &
    all(new_temperature == temperature), "invalid cut-state rollback")
  call advance_reactive_eb_redistributed_euler_3d( &
    species, state, temperature, geometry, "rusanov", 1.0e-8_dp, &
    new_state, new_temperature, ok)
  call require(.not. ok .and. all(new_state == state) .and. &
    all(new_temperature == temperature), &
    "redistributed invalid cut-state rollback")
  call advance_reactive_eb_state_redistributed_euler_3d( &
    species, state, temperature, geometry, "rusanov", 1.0e-8_dp, &
    new_state, new_temperature, ok)
  call require(.not. ok .and. all(new_state == state) .and. &
    all(new_temperature == temperature), &
    "StateRedist invalid cut-state rollback")

  write(*, '(a,es24.16)') &
    "maximum uniform state error: ", maximum_uniform_state_error
  write(*, '(a,es24.16)') &
    "maximum uniform temperature error: ", &
    maximum_uniform_temperature_error
  write(*, '(a,es24.16)') &
    "maximum conservation error: ", conservation_error
  write(*, '(a,es24.16)') &
    "maximum species closure error: ", species_closure_error
  write(*, '(a,es24.16)') &
    "small-cell full-grid timestep: ", small_cell_full_step
  write(*, '(a,es24.16)') &
    "FluxRedist cut-cell density: ", redistributed_cut_density
  write(*, '(a,es24.16)') &
    "StateRedist cut-cell density: ", state_redistributed_cut_density
  write(*, '(a)') "test_reactive_eb_hydro_3d: PASS"

contains

  subroutine verify_small_cell_hydro_relief()
    real(dp), allocatable :: raw_rhs(:, :, :, :)
    real(dp) :: stable_dt, local_conservation_error
    logical :: local_ok

    call build_axis_plane_eb_geometry_3d( &
      nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
      "x", 0.395_dp, geometry, local_ok)
    call require(local_ok, "small-cell hydro geometry")
    call assert_close(geometry%volume_fraction(4, 4, 3), 0.05_dp, &
      3.0e-13_dp, "hydro small-cell fraction")
    call build_cell_state(0.31_dp, [0.0_dp, 0.0_dp, 0.0_dp], &
      state_cell, base_temperature, sound_speed)
    call fill_state(state_cell, base_temperature)
    call build_cell_state(0.93_dp, [0.0_dp, 0.0_dp, 0.0_dp], &
      perturbed_cell, perturbed_temperature, perturbed_sound_speed)
    state(:, 4, 4, 3) = perturbed_cell
    temperature(4, 4, 3) = perturbed_temperature
    stable_dt = 0.5_dp / (max(sound_speed, perturbed_sound_speed) * &
      (1.0_dp / geometry%dx + 1.0_dp / geometry%dy + &
        1.0_dp / geometry%dz))
    small_cell_full_step = stable_dt

    call reactive_eb_outflow_riemann_fluxes_3d( &
      species, state, temperature, geometry, "rusanov", &
      x_flux, y_flux, z_flux, local_ok)
    call require(local_ok, "small-cell raw face fluxes")
    allocate(raw_rhs(nvar, nx, ny, nz))
    call reactive_eb_flux_divergence_3d( &
      species, state, temperature, geometry, x_flux, y_flux, z_flux, &
      raw_rhs, local_ok)
    call require(local_ok, "small-cell raw flux divergence")
    call require(state(irho, 4, 4, 3) + &
      stable_dt * raw_rhs(irho, 4, 4, 3) < 0.0_dp, &
      "raw full-cell timestep produces negative cut density")

    call advance_reactive_eb_euler_3d( &
      species, state, temperature, geometry, "rusanov", stable_dt, &
      new_state, new_temperature, local_ok)
    call require(.not. local_ok .and. all(new_state == state) .and. &
      all(new_temperature == temperature), &
      "raw small-cell hydro rejects full-cell timestep")

    call fluid_integrals(state, geometry, before_integrals)
    call advance_reactive_eb_redistributed_euler_3d( &
      species, state, temperature, geometry, "rusanov", stable_dt, &
      new_state, new_temperature, local_ok)
    call require(local_ok, &
      "redistribution accepts full-cell-scale hydro timestep")
    call require(new_state(irho, 4, 4, 3) > 0.0_dp, &
      "redistributed small-cell density remains positive")
    redistributed_cut_density = new_state(irho, 4, 4, 3)
    call fluid_integrals(new_state, geometry, after_integrals)
    local_conservation_error = maxval( &
      abs(after_integrals - before_integrals) / &
        max(1.0_dp, abs(before_integrals)))
    call require(local_conservation_error <= 3.0e-13_dp, &
      "small-cell full-step conservation")

    call advance_reactive_eb_state_redistributed_euler_3d( &
      species, state, temperature, geometry, "rusanov", stable_dt, &
      new_state, new_temperature, local_ok)
    call require(local_ok, &
      "StateRedist accepts full-cell-scale hydro timestep")
    call require(new_state(irho, 4, 4, 3) > 0.0_dp, &
      "StateRedist small-cell density remains positive")
    state_redistributed_cut_density = new_state(irho, 4, 4, 3)
    call fluid_integrals(new_state, geometry, after_integrals)
    local_conservation_error = maxval( &
      abs(after_integrals - before_integrals) / &
        max(1.0_dp, abs(before_integrals)))
    call require(local_conservation_error <= 3.0e-13_dp, &
      "StateRedist small-cell full-step conservation")
  end subroutine verify_small_cell_hydro_relief

  subroutine verify_uniform_tangent_update(axis, solver, velocity)
    character(len=*), intent(in) :: axis, solver
    real(dp), intent(in) :: velocity(3)

    real(dp) :: local_temperature, local_sound_speed
    real(dp) :: local_scale, local_error
    logical :: local_ok

    call build_axis_plane_eb_geometry_3d( &
      nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
      axis, 0.37_dp, geometry, local_ok)
    call require(local_ok, trim(axis)//"-plane uniform geometry")
    call build_cell_state( &
      0.31_dp, velocity, state_cell, local_temperature, local_sound_speed)
    call fill_state(state_cell, local_temperature)
    call advance_reactive_eb_euler_3d( &
      species, state, temperature, geometry, solver, 1.0e-7_dp, &
      new_state, new_temperature, local_ok)
    call require(local_ok, trim(axis)//"-plane uniform Euler update")
    local_scale = max(1.0_dp, maxval(abs(state)))
    local_error = maxval(abs(new_state - state)) / local_scale
    maximum_uniform_state_error = max( &
      maximum_uniform_state_error, local_error)
    call require(local_error <= 2.0e-14_dp, &
      trim(axis)//"-plane uniform state invariance")
    local_error = maxval(abs(new_temperature - temperature)) / &
      max(1.0_dp, local_temperature)
    maximum_uniform_temperature_error = max( &
      maximum_uniform_temperature_error, local_error)
    call require(local_error <= 2.0e-11_dp, &
      trim(axis)//"-plane uniform temperature invariance")
    call advance_reactive_eb_redistributed_euler_3d( &
      species, state, temperature, geometry, solver, 1.0e-7_dp, &
      new_state, new_temperature, local_ok)
    call require(local_ok, &
      trim(axis)//"-plane uniform redistributed Euler update")
    local_error = maxval(abs(new_state - state)) / local_scale
    maximum_uniform_state_error = max( &
      maximum_uniform_state_error, local_error)
    call require(local_error <= 2.0e-14_dp, &
      trim(axis)//"-plane redistributed state invariance")
    local_error = maxval(abs(new_temperature - temperature)) / &
      max(1.0_dp, local_temperature)
    maximum_uniform_temperature_error = max( &
      maximum_uniform_temperature_error, local_error)
    call require(local_error <= 2.0e-11_dp, &
      trim(axis)//"-plane redistributed temperature invariance")
    call advance_reactive_eb_state_redistributed_euler_3d( &
      species, state, temperature, geometry, solver, 1.0e-7_dp, &
      new_state, new_temperature, local_ok)
    call require(local_ok, &
      trim(axis)//"-plane uniform StateRedist Euler update")
    local_error = maxval(abs(new_state - state)) / local_scale
    maximum_uniform_state_error = max( &
      maximum_uniform_state_error, local_error)
    call require(local_error <= 2.0e-14_dp, &
      trim(axis)//"-plane StateRedist state invariance")
    local_error = maxval(abs(new_temperature - temperature)) / &
      max(1.0_dp, local_temperature)
    maximum_uniform_temperature_error = max( &
      maximum_uniform_temperature_error, local_error)
    call require(local_error <= 2.0e-11_dp, &
      trim(axis)//"-plane StateRedist temperature invariance")
  end subroutine verify_uniform_tangent_update

  subroutine build_cell_state( &
      density, velocity, conserved, recovered_temperature, recovered_sound)
    real(dp), intent(in) :: density, velocity(3)
    real(dp), intent(out) :: conserved(:)
    real(dp), intent(out) :: recovered_temperature, recovered_sound

    logical :: local_ok
    integer :: species_index

    primitive(1:5) = [density, velocity, pressure]
    do species_index = 1, size(species)
      primitive(reactive_mass_fraction_component(species_index)) = &
        mass_fractions(species_index)
    end do
    call reactive_primitive_to_conserved( &
      species, primitive, conserved, recovered_temperature, &
      recovered_sound, local_ok)
    call require(local_ok, "cell-state construction")
  end subroutine build_cell_state

  subroutine fill_state(conserved, cell_temperature)
    real(dp), intent(in) :: conserved(:), cell_temperature

    integer :: local_i, local_j, local_k

    do local_k = 1, nz
      do local_j = 1, ny
        do local_i = 1, nx
          state(:, local_i, local_j, local_k) = conserved
          temperature(local_i, local_j, local_k) = cell_temperature
        end do
      end do
    end do
  end subroutine fill_state

  subroutine fluid_integrals(field, local_geometry, integrals)
    real(dp), intent(in) :: field(:, :, :, :)
    type(eb_geometry_3d), intent(in) :: local_geometry
    real(dp), intent(out) :: integrals(:)

    real(dp) :: cell_volume
    integer :: component

    cell_volume = local_geometry%dx * local_geometry%dy * local_geometry%dz
    do component = 1, size(field, 1)
      integrals(component) = sum( &
        field(component, :, :, :) * local_geometry%volume_fraction) * &
        cell_volume
    end do
  end subroutine fluid_integrals

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

  subroutine assert_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: message

    if (abs(actual - expected) > tolerance) error stop message
  end subroutine assert_close

end program test_reactive_eb_hydro_3d
