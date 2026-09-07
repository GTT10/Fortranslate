module reactive_2d_application_mod
  use precision_mod, only: dp
  use constants_mod, only: pelef_version
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use simulation_config_reactive_2d_mod, only: reactive_2d_config
  use reactive_2d_mod, only: &
    simulate_reactive_2d, write_reactive_2d_csv, reactive_extrema_2d
  implicit none
  private

  public :: run_reactive_2d_application

contains

  subroutine run_reactive_2d_application( &
      input_path, application_label, bundle_sha256, config, species, &
      reactions, transport, base_mole_fractions, chemistry_integrator)
    character(len=*), intent(in) :: input_path, application_label
    character(len=*), intent(in) :: bundle_sha256
    type(reactive_2d_config), intent(in) :: config
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: base_mole_fractions(:)
    character(len=*), intent(in), optional :: chemistry_integrator

    real(dp), allocatable :: state(:, :, :), temperature(:, :)
    real(dp) :: dx, dy, time, base_density, minimum_theta
    real(dp) :: minimum_transport_theta
    real(dp) :: initial_integrals(5), final_integrals(5)
    real(dp) :: conservation_error(5)
    real(dp) :: minimum_density, maximum_density, minimum_pressure
    real(dp) :: maximum_pressure, minimum_temperature, maximum_temperature
    real(dp) :: maximum_speed, maximum_closure_error
    logical :: ok
    integer :: steps

    if (size(species) < 2 .or. size(reactions) < 1 .or. &
        size(transport) /= size(species) .or. &
        size(base_mole_fractions) /= size(species)) then
      error stop "Reactive 2D application data have incompatible dimensions"
    end if
    call simulate_reactive_2d( &
      species, reactions, config, state, temperature, dx, dy, time, steps, &
      initial_integrals, final_integrals, minimum_theta, base_density, ok, &
      transport, minimum_transport_theta, base_mole_fractions, &
      chemistry_integrator=chemistry_integrator)
    if (.not. ok) error stop "Reactive 2D simulation failed"
    call write_reactive_2d_csv( &
      config%output_file, species, config, state, temperature, dx, dy, time, ok)
    if (.not. ok) error stop "Reactive 2D output failed"
    call reactive_extrema_2d( &
      species, state, temperature, config%nx, config%ny, minimum_density, &
      maximum_density, minimum_pressure, maximum_pressure, &
      minimum_temperature, maximum_temperature, maximum_speed, &
      maximum_closure_error, ok)
    if (.not. ok) error stop "Reactive 2D diagnostics failed"

    conservation_error = abs(final_integrals - initial_integrals) / &
      max(1.0_dp, abs(initial_integrals))
    write(*, '(a)') &
      "PeleF " // pelef_version // " " // trim(application_label)
    if (len_trim(bundle_sha256) > 0) then
      write(*, '(a,1x,a)') "Bundle SHA-256:", trim(bundle_sha256)
      if (present(chemistry_integrator)) then
        write(*, '(a,1x,a)') "Chemistry integrator:", &
          trim(chemistry_integrator)
      end if
    end if
    write(*, '(a,1x,a)') "Input:", trim(input_path)
    write(*, '(a,1x,a)') "Problem:", trim(config%problem)
    write(*, '(a,i0,a,i0)') "Grid: ", config%nx, " x ", config%ny
    write(*, '(a,1x,a)') "Reconstruction:", trim(config%reconstruction)
    write(*, '(a,1x,a)') "Riemann solver:", trim(config%riemann_solver)
    write(*, '(a,l2)') "Transverse correction: ", &
      config%use_transverse_correction
    write(*, '(a,1x,a,1x,a)') "Boundary x:", &
      trim(config%boundary_x_lower), trim(config%boundary_x_upper)
    write(*, '(a,1x,a,1x,a)') "Boundary y:", &
      trim(config%boundary_y_lower), trim(config%boundary_y_upper)
    if (trim(config%reconstruction) == "characteristic_ppm") then
      write(*, '(a,l2)') "PPM contact steepening: ", &
        config%ppm_contact_steepening
      write(*, '(a,l2)') "PPM shock flattening: ", &
        config%ppm_shock_flattening
    end if
    write(*, '(a,l2)') "Chemistry: ", config%chemistry_enabled
    write(*, '(a,1x,a)') "Chemistry model:", trim(config%chemistry_model)
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
  end subroutine run_reactive_2d_application

end module reactive_2d_application_mod
