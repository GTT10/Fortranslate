program test_pelec_riemann
  use precision_mod, only: dp
  use state_indices_mod, only: &
    ncons, nprim, irho, imx, imy, imz, iet
  use state_conversion_mod, only: primitive_to_conserved
  use riemann_rusanov_mod, only: euler_physical_flux_x
  use riemann_pelec_mod, only: pelec_riemann_flux_x
  use riemann_flux_mod, only: compute_riemann_flux_x
  implicit none

  real(dp), parameter :: gamma = 1.4_dp
  real(dp), parameter :: tolerance = 5.0e-13_dp
  real(dp) :: left_primitive(nprim), right_primitive(nprim)
  real(dp) :: left_state(ncons), right_state(ncons)
  real(dp) :: physical_flux(ncons), numerical_flux(ncons), selected_flux(ncons)
  real(dp) :: rho_interface, velocity_interface, pressure_interface
  logical :: ok

  left_primitive = [1.0_dp, 0.5_dp, 0.2_dp, -0.1_dp, 1.0_dp]
  call primitive_to_conserved(left_primitive, gamma, left_state, ok)
  call assert_true(ok, "equal-state conversion")
  call euler_physical_flux_x(left_state, gamma, physical_flux, ok)
  call assert_true(ok, "physical flux")
  call pelec_riemann_flux_x( &
    left_state, left_state, gamma, numerical_flux, ok, &
    rho_interface, velocity_interface, pressure_interface)
  call assert_true(ok, "equal-state PeleC flux")
  call assert_close(maxval(abs(numerical_flux - physical_flux)), &
    0.0_dp, tolerance, "equal-state consistency")
  call assert_close(rho_interface, 1.0_dp, tolerance, "equal-state rho")
  call assert_close(velocity_interface, 0.5_dp, tolerance, "equal-state u")
  call assert_close(pressure_interface, 1.0_dp, tolerance, "equal-state p")

  left_primitive = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]
  right_primitive = [0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]
  call primitive_to_conserved(left_primitive, gamma, left_state, ok)
  call assert_true(ok, "contact left conversion")
  call primitive_to_conserved(right_primitive, gamma, right_state, ok)
  call assert_true(ok, "contact right conversion")
  call pelec_riemann_flux_x( &
    left_state, right_state, gamma, numerical_flux, ok, &
    rho_interface, velocity_interface, pressure_interface)
  call assert_true(ok, "stationary-contact flux")
  call assert_close(numerical_flux(irho), 0.0_dp, tolerance, "contact mass")
  call assert_close(numerical_flux(imx), 1.0_dp, tolerance, "contact momentum")
  call assert_close(numerical_flux(imy), 0.0_dp, tolerance, "contact y momentum")
  call assert_close(numerical_flux(imz), 0.0_dp, tolerance, "contact z momentum")
  call assert_close(numerical_flux(iet), 0.0_dp, tolerance, "contact energy")
  call assert_close(rho_interface, 0.75_dp, tolerance, "contact rho")
  call assert_close(velocity_interface, 0.0_dp, tolerance, "contact u")
  call assert_close(pressure_interface, 1.0_dp, tolerance, "contact p")

  left_primitive = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]
  right_primitive = [0.125_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.1_dp]
  call primitive_to_conserved(left_primitive, gamma, left_state, ok)
  call assert_true(ok, "Sod left conversion")
  call primitive_to_conserved(right_primitive, gamma, right_state, ok)
  call assert_true(ok, "Sod right conversion")
  call pelec_riemann_flux_x( &
    left_state, right_state, gamma, numerical_flux, ok, &
    rho_interface, velocity_interface, pressure_interface)
  call assert_true(ok, "Sod PeleC flux")
  call assert_close(rho_interface, 0.42178883109402565_dp, 2.0e-13_dp, &
    "Sod interface rho")
  call assert_close(velocity_interface, 0.6841486813454064_dp, 2.0e-13_dp, &
    "Sod interface u")
  call assert_close(pressure_interface, 0.19050436353163594_dp, 2.0e-13_dp, &
    "Sod interface p")
  call assert_close(numerical_flux(irho), 0.28856627259919798_dp, 5.0e-13_dp, &
    "Sod mass flux")
  call assert_close(numerical_flux(imx), 0.38792659841113630_dp, 5.0e-13_dp, &
    "Sod momentum flux")
  call assert_close(numerical_flux(iet), 0.52369966268303800_dp, 5.0e-13_dp, &
    "Sod energy flux")

  call compute_riemann_flux_x( &
    left_state, right_state, gamma, "pelec", selected_flux, ok)
  call assert_true(ok, "PeleC solver selection")
  call assert_close(maxval(abs(selected_flux - numerical_flux)), &
    0.0_dp, tolerance, "PeleC dispatch consistency")

  call compute_riemann_flux_x( &
    left_state, right_state, gamma, "unknown", selected_flux, ok)
  call assert_true(.not. ok, "unknown solver rejection")

  write(*, '(a)') "test_pelec_riemann: PASS"

contains

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*, '(a,1x,a)') "FAIL:", trim(label)
      error stop 1
    end if
  end subroutine assert_true

  subroutine assert_close(actual, expected, tol, label)
    real(dp), intent(in) :: actual, expected, tol
    character(len=*), intent(in) :: label
    if (abs(actual - expected) > tol) then
      write(*, '(a,1x,a,3(1x,es24.16))') &
        "FAIL:", trim(label), actual, expected, tol
      error stop 1
    end if
  end subroutine assert_close

end program test_pelec_riemann
