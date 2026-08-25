program pelef_reactive_eb_2d
  use precision_mod, only: dp
  use constants_mod, only: pelef_version
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
  use eb_geometry_2d_mod, only: &
    eb_geometry_2d, eb_covered_cell, eb_cut_cell, eb_regular_cell
  use simulation_config_reactive_eb_2d_mod, only: &
    reactive_eb_2d_config, read_reactive_eb_2d_configuration
  use reactive_eb_2d_driver_mod, only: &
    simulate_reactive_eb_2d, write_reactive_eb_2d_csv, &
    reactive_eb_extrema_2d
  implicit none

  type(reactive_eb_2d_config) :: config
  type(eb_geometry_2d) :: geometry
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  real(dp), allocatable :: state(:, :, :), temperature(:, :)
  real(dp), allocatable :: initial_integrals(:), final_integrals(:)
  real(dp) :: time, minimum_dt, base_density, conservation_error
  real(dp) :: minimum_transport_theta
  real(dp) :: minimum_density, maximum_density
  real(dp) :: minimum_pressure, maximum_pressure
  real(dp) :: minimum_temperature, maximum_temperature
  real(dp) :: maximum_speed, maximum_closure_error
  character(len=1024) :: input_path, message
  logical :: ok
  integer :: steps

  if (command_argument_count() /= 1) then
    write(*, '(a)') "Usage: pelef_reactive_eb_2d <input.nml>"
    error stop 2
  end if
  call get_command_argument(1, input_path)
  call read_reactive_eb_2d_configuration( &
    trim(input_path), config, ok, message)
  if (.not. ok) then
    write(*, '(a)') trim(message)
    error stop 2
  end if

  select case (trim(config%flow%chemistry_model))
  case ("elementary")
    call load_h2o2_elementary_thermo(species, ok)
    if (.not. ok) error stop "Failed to load elementary thermodynamics"
    call load_h2o2_elementary_mechanism(reactions, ok)
    if (.not. ok) error stop "Failed to load elementary mechanism"
    call load_h2o2_elementary_transport(transport, ok)
    if (.not. ok) error stop "Failed to load elementary transport"
  case ("full_h2o2")
    call load_h2o2_full_thermo(species, ok)
    if (.not. ok) error stop "Failed to load full H2/O2 thermodynamics"
    call load_h2o2_full_mechanism(reactions, ok)
    if (.not. ok) error stop "Failed to load full H2/O2 mechanism"
    call load_h2o2_full_transport(transport, ok)
    if (.not. ok) error stop "Failed to load full H2/O2 transport"
  case default
    error stop "Unknown chemistry model"
  end select

  call simulate_reactive_eb_2d( &
    species, reactions, config, state, temperature, geometry, time, steps, &
    initial_integrals, final_integrals, minimum_dt, base_density, ok, &
    transport, minimum_transport_theta)
  if (.not. ok) error stop "Reactive EB 2D simulation failed"
  call write_reactive_eb_2d_csv( &
    config%flow%output_file, species, config, state, temperature, &
    geometry, time, ok)
  if (.not. ok) error stop "Reactive EB 2D output failed"
  call reactive_eb_extrema_2d( &
    species, state, temperature, geometry, minimum_density, &
    maximum_density, minimum_pressure, maximum_pressure, &
    minimum_temperature, maximum_temperature, maximum_speed, &
    maximum_closure_error, ok)
  if (.not. ok) error stop "Reactive EB 2D diagnostics failed"

  conservation_error = maxval(abs(final_integrals - initial_integrals) / &
    max(1.0_dp, abs(initial_integrals)))
  write(*, '(a)') "PeleF " // pelef_version // " reactive EB 2D"
  write(*, '(a,1x,a)') "Problem:", trim(config%flow%problem)
  write(*, '(a,1x,a)') "Geometry:", trim(config%geometry)
  write(*, '(a,i0,a,i0)') "Grid: ", geometry%nx, " x ", geometry%ny
  write(*, '(a,1x,a)') "Riemann solver:", &
    trim(config%flow%riemann_solver)
  write(*, '(a,1x,a)') "Reconstruction:", &
    trim(config%flow%reconstruction)
  write(*, '(a,1x,a)') "Limiter:", trim(config%flow%limiter)
  write(*, '(a,l2)') "Chemistry: ", config%flow%chemistry_enabled
  write(*, '(a,1x,a)') "Chemistry model:", &
    trim(config%flow%chemistry_model)
  write(*, '(a,l2)') "Molecular transport: ", &
    config%flow%transport_enabled
  if (config%flow%transport_enabled) then
    write(*, '(a,l2)') "Viscosity: ", config%flow%viscosity_enabled
    write(*, '(a,l2)') "Thermal conduction: ", &
      config%flow%thermal_conduction_enabled
    write(*, '(a,l2)') "Species diffusion: ", &
      config%flow%species_diffusion_enabled
    write(*, '(a,l2)') "Barodiffusion: ", &
      config%flow%barodiffusion_enabled
    write(*, '(a,es24.16)') "Minimum transport limiter theta: ", &
      minimum_transport_theta
  end if
  write(*, '(a,es24.16)') "StateRedist target volume fraction: ", &
    config%state_redist_target_volume_fraction
  write(*, '(a,i0)') "StateRedist max order: ", &
    config%state_redist_max_order
  write(*, '(a,i0)') "Regular cells: ", &
    count(geometry%cell_type == eb_regular_cell)
  write(*, '(a,i0)') "Cut cells: ", &
    count(geometry%cell_type == eb_cut_cell)
  write(*, '(a,i0)') "Covered cells: ", &
    count(geometry%cell_type == eb_covered_cell)
  write(*, '(a,i0)') "Completed steps: ", steps
  write(*, '(a,es24.16)') "Final time: ", time
  write(*, '(a,es24.16)') "Minimum accepted dt: ", minimum_dt
  write(*, '(a,es24.16)') "Maximum conservation error: ", &
    conservation_error
  write(*, '(a,es24.16)') "Minimum density: ", minimum_density
  write(*, '(a,es24.16)') "Minimum pressure: ", minimum_pressure
  write(*, '(a,es24.16)') "Temperature range: ", &
    maximum_temperature - minimum_temperature
  write(*, '(a,es24.16)') "Maximum speed: ", maximum_speed
  write(*, '(a,es24.16)') "Maximum composition closure error: ", &
    maximum_closure_error
  write(*, '(a,1x,a)') "Output:", trim(config%flow%output_file)
end program pelef_reactive_eb_2d
