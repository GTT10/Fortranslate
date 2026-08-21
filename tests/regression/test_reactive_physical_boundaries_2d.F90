program test_reactive_physical_boundaries_2d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use transport_database_mod, only: &
    gas_transport_species, load_h2o2_elementary_transport
  use simulation_config_reactive_2d_mod, only: reactive_2d_config
  use reactive_1d_mod, only: reactive_nprim, reactive_conserved_to_primitive
  use reactive_2d_mod, only: simulate_reactive_2d
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  logical :: ok

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "thermodynamics load")
  call load_h2o2_elementary_mechanism(reactions, ok)
  call require(ok, "mechanism load")
  call load_h2o2_elementary_transport(transport, ok)
  call require(ok, "transport load")
  call test_couette_profile()
  call test_thermal_walls()
  call test_uniform_inflow_outflow()

contains

  subroutine base_config(config)
    type(reactive_2d_config), intent(out) :: config
    config%nx = 8
    config%ny = 16
    config%x_lower = 0.0_dp
    config%x_upper = 0.004_dp
    config%y_lower = 0.0_dp
    config%y_upper = 0.004_dp
    config%final_time = 1.0e-6_dp
    config%cfl = 0.20_dp
    config%maximum_steps = 20000
    config%reconstruction = "pcm"
    config%riemann_solver = "hllc"
    config%limiter = "mc"
    config%use_transverse_correction = .false.
    config%chemistry_enabled = .false.
    config%transport_enabled = .true.
    config%transport_cfl = 0.35_dp
    config%initial_temperature = 1000.0_dp
    config%initial_pressure = 101325.0_dp
    config%initial_velocity_x = 0.0_dp
    config%initial_velocity_y = 0.0_dp
    config%x_h2 = 0.29570_dp
    config%x_h = 1.0e-5_dp
    config%x_o = 1.0e-5_dp
    config%x_o2 = 0.14784_dp
    config%x_oh = 1.0e-5_dp
    config%x_h2o = 0.0_dp
    config%x_n2 = 0.55643_dp
  end subroutine base_config

  subroutine run_case(config, state, temperature, dx, dy)
    type(reactive_2d_config), intent(in) :: config
    real(dp), allocatable, intent(out) :: state(:, :, :), temperature(:, :)
    real(dp), intent(out) :: dx, dy
    real(dp) :: time, initial_integrals(5), final_integrals(5)
    real(dp) :: theta, transport_theta, base_density
    logical :: local_ok
    integer :: steps
    call simulate_reactive_2d( &
      species, reactions, config, state, temperature, dx, dy, time, steps, &
      initial_integrals, final_integrals, theta, base_density, local_ok, &
      transport, transport_theta)
    call require(local_ok, "physical-boundary simulation")
    call require(steps > 0, "physical-boundary simulation advances")
  end subroutine run_case

  subroutine test_couette_profile()
    type(reactive_2d_config) :: config
    real(dp), allocatable :: state(:, :, :), temperature(:, :), primitive(:)
    real(dp) :: dx, dy, expected, local_t, sound, error
    logical :: local_ok
    integer :: i, j

    call base_config(config)
    config%problem = "couette_channel"
    config%viscosity_enabled = .true.
    config%thermal_conduction_enabled = .false.
    config%species_diffusion_enabled = .false.
    config%barodiffusion_enabled = .false.
    config%boundary_y_lower = "no_slip_wall"
    config%boundary_y_upper = "no_slip_wall"
    config%wall_velocity_y_lower = [0.0_dp, 0.0_dp, 0.0_dp]
    config%wall_velocity_y_upper = [20.0_dp, 0.0_dp, 0.0_dp]
    call run_case(config, state, temperature, dx, dy)
    allocate(primitive(reactive_nprim(size(species))))
    error = 0.0_dp
    do j = 1, config%ny
      expected = 20.0_dp * (real(j, dp) - 0.5_dp) / real(config%ny, dp)
      do i = 1, config%nx
        call reactive_conserved_to_primitive( &
          species, state(:, i, j), temperature(i, j), primitive, local_t, &
          sound, local_ok)
        call require(local_ok, "Couette primitive recovery")
        error = max(error, abs(primitive(2) - expected))
        call require(abs(primitive(3)) < 2.0e-9_dp, "Couette no penetration")
      end do
    end do
    call require(error < 2.0e-8_dp, "linear Couette profile is preserved")
  end subroutine test_couette_profile

  subroutine test_thermal_walls()
    type(reactive_2d_config) :: config
    real(dp), allocatable :: state(:, :, :), temperature(:, :)
    real(dp) :: dx, dy
    integer :: j

    call base_config(config)
    config%problem = "thermal_channel"
    config%final_time = 5.0e-7_dp
    config%viscosity_enabled = .false.
    config%thermal_conduction_enabled = .true.
    config%species_diffusion_enabled = .false.
    config%barodiffusion_enabled = .false.
    config%boundary_y_lower = "slip_wall"
    config%boundary_y_upper = "slip_wall"
    config%thermal_y_lower = "isothermal"
    config%thermal_y_upper = "isothermal"
    config%wall_temperature_y_lower = 800.0_dp
    config%wall_temperature_y_upper = 1200.0_dp
    call run_case(config, state, temperature, dx, dy)
    call require(minval(temperature) > 790.0_dp, "thermal wall lower bound")
    call require(maxval(temperature) < 1210.0_dp, "thermal wall upper bound")
    do j = 1, config%ny - 1
      call require(sum(temperature(:, j + 1)) >= sum(temperature(:, j)), &
        "thermal channel remains monotone")
    end do
  end subroutine test_thermal_walls

  subroutine test_uniform_inflow_outflow()
    type(reactive_2d_config) :: config
    real(dp), allocatable :: state(:, :, :), temperature(:, :), primitive(:)
    real(dp) :: dx, dy, local_t, sound, maximum_error
    logical :: local_ok
    integer :: i, j

    call base_config(config)
    config%nx = 16
    config%ny = 6
    config%problem = "inflow_outflow"
    config%reconstruction = "characteristic_plm"
    config%use_transverse_correction = .true.
    config%initial_velocity_x = 75.0_dp
    config%boundary_x_lower = "inflow"
    config%boundary_x_upper = "outflow"
    config%viscosity_enabled = .true.
    config%thermal_conduction_enabled = .true.
    config%species_diffusion_enabled = .true.
    config%barodiffusion_enabled = .true.
    call run_case(config, state, temperature, dx, dy)
    allocate(primitive(reactive_nprim(size(species))))
    maximum_error = 0.0_dp
    do j = 1, config%ny
      do i = 1, config%nx
        call reactive_conserved_to_primitive( &
          species, state(:, i, j), temperature(i, j), primitive, local_t, &
          sound, local_ok)
        call require(local_ok, "inflow primitive recovery")
        maximum_error = max(maximum_error, abs(primitive(2) - 75.0_dp))
        maximum_error = max(maximum_error, abs(primitive(3)))
        maximum_error = max(maximum_error, abs(primitive(5) - 101325.0_dp) / 101325.0_dp)
        maximum_error = max(maximum_error, abs(local_t - 1000.0_dp) / 1000.0_dp)
      end do
    end do
    call require(maximum_error < 5.0e-12_dp, &
      "uniform fixed inflow/outflow state is preserved")
  end subroutine test_uniform_inflow_outflow

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') "FAILED: " // trim(message)
      error stop 1
    end if
  end subroutine require
end program test_reactive_physical_boundaries_2d
