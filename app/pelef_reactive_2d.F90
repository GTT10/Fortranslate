program pelef_reactive_2d
  use precision_mod, only: dp
  use constants_mod, only: pelef_version
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use transport_database_mod, only: &
    gas_transport_species, load_h2o2_elementary_transport
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use simulation_config_reactive_2d_mod, only: &
    reactive_2d_config, read_reactive_2d_configuration
  use reactive_2d_mod, only: &
    simulate_reactive_2d, write_reactive_2d_csv, reactive_extrema_2d
  implicit none

  type(reactive_2d_config) :: config
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  real(dp), allocatable :: state(:, :, :), temperature(:, :)
  real(dp) :: dx, dy, time, base_density, minimum_theta
  real(dp) :: minimum_transport_theta
  real(dp) :: initial_integrals(5), final_integrals(5), conservation_error(5)
  real(dp) :: minimum_density, maximum_density, minimum_pressure
  real(dp) :: maximum_pressure, minimum_temperature, maximum_temperature
  real(dp) :: maximum_speed, maximum_closure_error
  character(len=1024) :: input_path, message
  logical :: ok
  integer :: steps

  if (command_argument_count() /= 1) then
    write(*, '(a)') "Usage: pelef_reactive_2d <input.nml>"
    error stop 2
  end if
  call get_command_argument(1, input_path)
  call read_reactive_2d_configuration(trim(input_path), config, ok, message)
  if (.not. ok) then
    write(*, '(a)') trim(message)
    error stop 2
  end if
  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok) error stop "Failed to load reactive thermodynamics"
  call load_h2o2_elementary_mechanism(reactions, ok)
  if (.not. ok) error stop "Failed to load reactive mechanism"
  call load_h2o2_elementary_transport(transport, ok)
  if (.not. ok) error stop "Failed to load reactive transport database"

  call simulate_reactive_2d( &
    species, reactions, config, state, temperature, dx, dy, time, steps, &
    initial_integrals, final_integrals, minimum_theta, base_density, ok, &
    transport, minimum_transport_theta)
  if (.not. ok) error stop "Reactive 2D simulation failed"
  call write_reactive_2d_csv( &
    config%output_file, species, config, state, temperature, dx, dy, time, ok)
  if (.not. ok) error stop "Reactive 2D output failed"
  call reactive_extrema_2d( &
    species, state, temperature, config%nx, config%ny, minimum_density, &
    maximum_density, minimum_pressure, maximum_pressure, minimum_temperature, &
    maximum_temperature, maximum_speed, maximum_closure_error, ok)
  if (.not. ok) error stop "Reactive 2D diagnostics failed"

  conservation_error = abs(final_integrals - initial_integrals) / &
    max(1.0_dp, abs(initial_integrals))
  write(*, '(a)') "PeleF " // pelef_version // " reactive 2D"
  write(*, '(a,1x,a)') "Problem:", trim(config%problem)
  write(*, '(a,i0,a,i0)') "Grid: ", config%nx, " x ", config%ny
  write(*, '(a,1x,a)') "Reconstruction:", trim(config%reconstruction)
  write(*, '(a,1x,a)') "Riemann solver:", trim(config%riemann_solver)
  write(*, '(a,l2)') "Transverse correction: ", &
    config%use_transverse_correction
  if (trim(config%reconstruction) == "characteristic_ppm") then
    write(*, '(a,l2)') "PPM contact steepening: ", &
      config%ppm_contact_steepening
    write(*, '(a,l2)') "PPM shock flattening: ", &
      config%ppm_shock_flattening
  end if
  write(*, '(a,l2)') "Chemistry: ", config%chemistry_enabled
  write(*, '(a,l2)') "Molecular transport: ", config%transport_enabled
  if (config%transport_enabled) then
    write(*, '(a,l2)') "Viscosity: ", config%viscosity_enabled
    write(*, '(a,l2)') "Thermal conduction: ", &
      config%thermal_conduction_enabled
    write(*, '(a,l2)') "Species diffusion: ", &
      config%species_diffusion_enabled
    write(*, '(a,l2)') "Barodiffusion: ", config%barodiffusion_enabled
  end if
  write(*, '(a,i0)') "Completed steps: ", steps
  write(*, '(a,es24.16)') "Final time: ", time
  write(*, '(a,es24.16)') "Minimum transverse theta: ", minimum_theta
  write(*, '(a,es24.16)') "Minimum transport theta: ", &
    minimum_transport_theta
  write(*, '(a,es24.16)') "Maximum conservation error: ", &
    maxval(conservation_error)
  write(*, '(a,es24.16)') "Minimum density: ", minimum_density
  write(*, '(a,es24.16)') "Minimum pressure: ", minimum_pressure
  write(*, '(a,es24.16)') "Temperature range: ", &
    maximum_temperature - minimum_temperature
  write(*, '(a,es24.16)') "Maximum speed: ", maximum_speed
  write(*, '(a,es24.16)') "Maximum composition closure error: ", &
    maximum_closure_error
  write(*, '(a,1x,a)') "Output:", trim(config%output_file)
end program pelef_reactive_2d
