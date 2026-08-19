program test_reactive_hllc
  use precision_mod, only: dp
  use state_indices_mod, only: irho, imx, imy, imz, iet
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: &
    mass_fractions_from_mole_fractions, mixture_density
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_species_component, reactive_primitive_to_conserved, &
    reactive_hllc_flux_x, reactive_riemann_flux_x
  implicit none

  type(nasa7_species), allocatable :: species(:)
  real(dp), allocatable :: ql(:), qr(:), ul(:), ur(:), flux(:), expected(:)
  real(dp) :: xl(7), xr(7), yl(7), yr(7)
  real(dp) :: pressure, temperature, velocity, rho_l, rho_r
  real(dp) :: tl, tr, cl, cr
  logical :: ok
  integer :: k

  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok) error stop "Failed to load HLLC thermodynamics"
  allocate(ql(reactive_nprim(7)), qr(reactive_nprim(7)))
  allocate(ul(reactive_nvar(7)), ur(reactive_nvar(7)))
  allocate(flux(reactive_nvar(7)), expected(reactive_nvar(7)))

  xl = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
    1.0e-5_dp, 0.0_dp, 0.55643_dp]
  xr = xl
  xr(1) = xr(1) + 0.05_dp
  xr(7) = xr(7) - 0.05_dp
  call mass_fractions_from_mole_fractions(species, xl, yl, ok)
  if (.not. ok) error stop "Failed to construct left HLLC composition"
  call mass_fractions_from_mole_fractions(species, xr, yr, ok)
  if (.not. ok) error stop "Failed to construct right HLLC composition"

  pressure = 1.5e5_dp
  temperature = 1000.0_dp
  rho_l = mixture_density(species, yl, pressure, temperature, ok)
  if (.not. ok) error stop "Failed to construct left HLLC density"
  rho_r = mixture_density(species, yr, pressure, temperature, ok)
  if (.not. ok) error stop "Failed to construct right HLLC density"

  velocity = 0.0_dp
  call make_state(rho_l, velocity, pressure, yl, ql, ul, tl, cl)
  call make_state(rho_r, velocity, pressure, yr, qr, ur, tr, cr)
  call reactive_hllc_flux_x(species, ul, ur, tl, tr, flux, ok)
  if (.not. ok) error stop "Stationary material-contact HLLC flux failed"
  expected = 0.0_dp
  expected(imx) = pressure
  call assert_vector_close(flux, expected, 1.0e-7_dp, &
    "stationary material contact")

  velocity = 80.0_dp
  call make_state(rho_l, velocity, pressure, yl, ql, ul, tl, cl)
  call make_state(rho_r, velocity, pressure, yr, qr, ur, tr, cr)
  call reactive_hllc_flux_x(species, ul, ur, tl, tr, flux, ok)
  if (.not. ok) error stop "Moving material-contact HLLC flux failed"
  expected = 0.0_dp
  expected(irho) = rho_l * velocity
  expected(imx) = rho_l * velocity**2 + pressure
  expected(imy) = 0.0_dp
  expected(imz) = 0.0_dp
  expected(iet) = (ul(iet) + pressure) * velocity
  do k = 1, 7
    expected(reactive_species_component(k)) = rho_l * velocity * yl(k)
  end do
  expected(reactive_species_component(7)) = expected(irho) - &
    sum(expected(reactive_species_component(1): &
      reactive_species_component(6)))
  call assert_vector_close(flux, expected, 2.0e-10_dp, &
    "moving material contact")
  call assert_close(sum(flux(6:12)), flux(irho), 2.0e-12_dp, &
    "HLLC species closure")

  call reactive_riemann_flux_x(species, ul, ur, tl, tr, "hllc", flux, ok)
  if (.not. ok) error stop "HLLC dispatcher failed"
  call reactive_riemann_flux_x(species, ul, ur, tl, tr, "unknown", flux, ok)
  if (ok) error stop "Unknown reactive Riemann solver was accepted"

  write(*, '(a)') "test_reactive_hllc: PASS"

contains

  subroutine make_state(rho, u, p, y, q, conserved, t, c)
    real(dp), intent(in) :: rho, u, p, y(:)
    real(dp), intent(out) :: q(:), conserved(:), t, c
    logical :: local_ok
    integer :: species_index

    q(1:5) = [rho, u, 0.0_dp, 0.0_dp, p]
    do species_index = 1, size(y)
      q(reactive_mass_fraction_component(species_index)) = y(species_index)
    end do
    call reactive_primitive_to_conserved( &
      species, q, conserved, t, c, local_ok)
    if (.not. local_ok) error stop "Failed to construct HLLC state"
  end subroutine make_state

  subroutine assert_close(actual, expected_value, relative_tolerance, label)
    real(dp), intent(in) :: actual, expected_value, relative_tolerance
    character(len=*), intent(in) :: label
    real(dp) :: error

    error = abs(actual - expected_value) / max(1.0_dp, abs(expected_value))
    if (error > relative_tolerance) then
      write(*, '(a,3(1x,es24.16))') &
        trim(label), actual, expected_value, error
      error stop "Reactive HLLC scalar mismatch"
    end if
  end subroutine assert_close

  subroutine assert_vector_close(actual, expected_value, relative_tolerance, label)
    real(dp), intent(in) :: actual(:), expected_value(:), relative_tolerance
    character(len=*), intent(in) :: label
    real(dp) :: error

    error = maxval(abs(actual - expected_value) / max(1.0_dp, abs(expected_value)))
    if (error > relative_tolerance) then
      write(*, '(a,1x,es24.16)') trim(label), error
      error stop "Reactive HLLC vector mismatch"
    end if
  end subroutine assert_vector_close

end program test_reactive_hllc
