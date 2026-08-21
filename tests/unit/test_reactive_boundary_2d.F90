program test_reactive_boundary_2d
  use precision_mod, only: dp
  use state_indices_mod, only: imx, iet
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions, mixture_density
  use transport_database_mod, only: &
    gas_transport_species, load_h2o2_elementary_transport
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_species_component, &
    reactive_mass_fraction_component, reactive_primitive_to_conserved
  use reactive_boundary_2d_mod, only: &
    reactive_boundary_set_2d, initialize_periodic_boundary_set_2d, &
    sample_reactive_primitive_2d, boundary_x_lower, boundary_x_upper, &
    boundary_y_lower, boundary_y_upper
  use reactive_transport_2d_mod, only: reactive_transport_fluxes_2d_faces
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  logical :: ok

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "thermodynamics load")
  call load_h2o2_elementary_transport(transport, ok)
  call require(ok, "transport database load")
  call test_ghost_sampling()
  call test_wall_transport_fluxes()

contains

  subroutine make_uniform_state(nx, ny, primitive, state, temperature)
    integer, intent(in) :: nx, ny
    real(dp), intent(out) :: primitive(:)
    real(dp), intent(out) :: state(:, :, :), temperature(:, :)
    real(dp) :: xmol(7), y(7), rho, sound
    logical :: local_ok
    integer :: i, j, k

    xmol = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
      1.0e-5_dp, 0.0_dp, 0.55643_dp]
    call mass_fractions_from_mole_fractions(species, xmol, y, local_ok)
    call require(local_ok, "composition")
    rho = mixture_density(species, y, 101325.0_dp, 1000.0_dp, local_ok)
    call require(local_ok, "density")
    primitive(1:5) = [rho, 3.0_dp, 0.0_dp, 2.0_dp, 101325.0_dp]
    do k = 1, size(species)
      primitive(reactive_mass_fraction_component(k)) = y(k)
    end do
    do j = 1, ny
      do i = 1, nx
        call reactive_primitive_to_conserved( &
          species, primitive, state(:, i, j), temperature(i, j), sound, local_ok)
        call require(local_ok, "state construction")
      end do
    end do
  end subroutine make_uniform_state

  subroutine test_ghost_sampling()
    integer, parameter :: nx = 4, ny = 4
    type(reactive_boundary_set_2d) :: boundaries
    real(dp), allocatable :: primitive(:, :, :), state(:, :, :), temperature(:, :)
    real(dp), allocatable :: base(:), sampled(:)
    real(dp) :: sampled_temperature
    logical :: local_ok

    allocate(base(reactive_nprim(size(species))))
    allocate(sampled(size(base)))
    allocate(primitive(size(base), nx, ny))
    allocate(state(reactive_nvar(size(species)), nx, ny), temperature(nx, ny))
    call make_uniform_state(nx, ny, base, state, temperature)
    primitive = spread(spread(base, 2, nx), 3, ny)
    call initialize_periodic_boundary_set_2d(size(base), boundaries)
    boundaries%face(boundary_y_lower)%kind = "no_slip_wall"
    boundaries%face(boundary_y_upper)%kind = "no_slip_wall"
    boundaries%face(boundary_y_lower)%thermal = "isothermal"
    boundaries%face(boundary_y_lower)%wall_temperature = 800.0_dp
    boundaries%face(boundary_y_lower)%wall_velocity = [10.0_dp, 0.0_dp, 0.0_dp]
    call sample_reactive_primitive_2d( &
      primitive, temperature, nx, ny, 2, 0, boundaries, sampled, &
      sampled_temperature, local_ok)
    call require(local_ok, "no-slip ghost sample")
    call require(abs(sampled(2) - 17.0_dp) < 1.0e-13_dp, &
      "no-slip tangential reflection")
    call require(abs(sampled(3)) < 1.0e-13_dp, "wall normal reflection")
    call require(abs(sampled(4) + 2.0_dp) < 1.0e-13_dp, &
      "out-of-plane no-slip reflection")
    call require(abs(sampled_temperature - 600.0_dp) < 1.0e-12_dp, &
      "isothermal temperature reflection")

    boundaries%face(boundary_y_lower)%kind = "slip_wall"
    call sample_reactive_primitive_2d( &
      primitive, temperature, nx, ny, 2, 0, boundaries, sampled, &
      sampled_temperature, local_ok)
    call require(local_ok, "slip ghost sample")
    call require(abs(sampled(2) - base(2)) < 1.0e-13_dp, &
      "slip tangential velocity")
    call require(abs(sampled(4) - base(4)) < 1.0e-13_dp, &
      "slip out-of-plane velocity")

    boundaries%face(boundary_x_lower)%kind = "inflow"
    boundaries%face(boundary_x_upper)%kind = "outflow"
    boundaries%face(boundary_x_lower)%inflow_primitive = base
    boundaries%face(boundary_x_lower)%inflow_primitive(2) = 55.0_dp
    boundaries%face(boundary_x_lower)%inflow_temperature = 900.0_dp
    call sample_reactive_primitive_2d( &
      primitive, temperature, nx, ny, 0, 2, boundaries, sampled, &
      sampled_temperature, local_ok)
    call require(local_ok, "inflow sample")
    call require(abs(sampled(2) - 55.0_dp) < 1.0e-13_dp, "fixed inflow state")
    call require(abs(sampled_temperature - 900.0_dp) < 1.0e-13_dp, &
      "fixed inflow temperature")
    call sample_reactive_primitive_2d( &
      primitive, temperature, nx, ny, nx + 1, 2, boundaries, sampled, &
      sampled_temperature, local_ok)
    call require(local_ok, "outflow sample")
    call require(maxval(abs(sampled - base)) < 1.0e-13_dp, &
      "outflow constant extrapolation")
  end subroutine test_ghost_sampling

  subroutine test_wall_transport_fluxes()
    integer, parameter :: nx = 4, ny = 4
    type(reactive_boundary_set_2d) :: boundaries
    real(dp), allocatable :: primitive(:), state(:, :, :), temperature(:, :)
    real(dp), allocatable :: flux_x(:, :, :), flux_y(:, :, :)
    real(dp) :: theta
    logical :: local_ok
    integer :: k

    allocate(primitive(reactive_nprim(size(species))))
    allocate(state(reactive_nvar(size(species)), nx, ny), temperature(nx, ny))
    allocate(flux_x(reactive_nvar(size(species)), 0:nx, ny))
    allocate(flux_y(reactive_nvar(size(species)), nx, 0:ny))
    call make_uniform_state(nx, ny, primitive, state, temperature)
    call initialize_periodic_boundary_set_2d(size(primitive), boundaries)
    boundaries%face(boundary_y_lower)%kind = "no_slip_wall"
    boundaries%face(boundary_y_upper)%kind = "no_slip_wall"
    boundaries%face(boundary_y_lower)%wall_velocity = [10.0_dp, 0.0_dp, 0.0_dp]
    boundaries%face(boundary_y_upper)%wall_velocity = [3.0_dp, 0.0_dp, 2.0_dp]
    call reactive_transport_fluxes_2d_faces( &
      species, transport, state, temperature, nx, ny, 1.0e-3_dp, 1.0e-3_dp, &
      1.0e-7_dp, .true., .false., .true., .true., flux_x, flux_y, theta, &
      local_ok, boundaries)
    call require(local_ok, "wall viscous/species flux")
    call require(flux_y(imx, 1, 0) > 0.0_dp, &
      "moving lower wall transfers positive x momentum")
    do k = 1, size(species)
      call require(abs(flux_y(reactive_species_component(k), 1, 0)) < &
        1.0e-20_dp, "impermeable wall species flux")
    end do

    boundaries%face(boundary_y_lower)%kind = "slip_wall"
    boundaries%face(boundary_y_upper)%kind = "slip_wall"
    boundaries%face(boundary_y_lower)%thermal = "adiabatic"
    boundaries%face(boundary_y_upper)%thermal = "adiabatic"
    call reactive_transport_fluxes_2d_faces( &
      species, transport, state, temperature, nx, ny, 1.0e-3_dp, 1.0e-3_dp, &
      1.0e-7_dp, .false., .true., .false., .false., flux_x, flux_y, theta, &
      local_ok, boundaries)
    call require(local_ok, "adiabatic wall flux")
    call require(abs(flux_y(iet, 1, 0)) < 1.0e-12_dp, &
      "adiabatic wall heat flux is zero")

    boundaries%face(boundary_y_lower)%thermal = "isothermal"
    boundaries%face(boundary_y_lower)%wall_temperature = 800.0_dp
    call reactive_transport_fluxes_2d_faces( &
      species, transport, state, temperature, nx, ny, 1.0e-3_dp, 1.0e-3_dp, &
      1.0e-7_dp, .false., .true., .false., .false., flux_x, flux_y, theta, &
      local_ok, boundaries)
    call require(local_ok, "isothermal wall flux")
    call require(flux_y(iet, 1, 0) < 0.0_dp, &
      "cold lower wall removes energy")
  end subroutine test_wall_transport_fluxes

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') "FAILED: " // trim(message)
      error stop 1
    end if
  end subroutine require
end program test_reactive_boundary_2d
