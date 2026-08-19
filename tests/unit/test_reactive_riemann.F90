program test_reactive_riemann
  use precision_mod, only: dp
  use state_indices_mod, only: irho, imx, imy, imz, iet
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: &
    mass_fractions_from_mole_fractions, mixture_density
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_primitive_to_conserved, reactive_physical_flux_x, &
    reactive_rusanov_flux_x, reactive_hllc_flux_x, &
    reactive_pelec_flux_x, &
    compute_reactive_riemann_flux_x
  implicit none

  type(nasa7_species), allocatable :: species(:)
  real(dp), allocatable :: ql(:), qr(:), ul(:), ur(:)
  real(dp), allocatable :: flux(:), pelec(:), reference(:), rusanov(:), primitive(:)
  real(dp) :: xl(7), xr(7), yl(7), yr(7)
  real(dp) :: rho_l, rho_r, temperature, pressure, velocity
  real(dp) :: tl, tr, cl, cr, dummy_t, dummy_c, error
  logical :: ok
  integer :: k

  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok) error stop "Failed to load HLLC thermodynamics"
  allocate(ql(reactive_nprim(7)), qr(reactive_nprim(7)))
  allocate(ul(reactive_nvar(7)), ur(reactive_nvar(7)))
  allocate(flux(reactive_nvar(7)), pelec(reactive_nvar(7)))
  allocate(reference(reactive_nvar(7)), rusanov(reactive_nvar(7)))
  allocate(primitive(reactive_nprim(7)))

  temperature = 1000.0_dp
  pressure = 101325.0_dp
  xl = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]
  xr = [0.25_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.75_dp]
  call mass_fractions_from_mole_fractions(species, xl, yl, ok)
  if (.not. ok) error stop "Failed to convert left HLLC composition"
  call mass_fractions_from_mole_fractions(species, xr, yr, ok)
  if (.not. ok) error stop "Failed to convert right HLLC composition"
  rho_l = mixture_density(species, yl, pressure, temperature, ok)
  if (.not. ok) error stop "Failed to evaluate left HLLC density"
  rho_r = mixture_density(species, yr, pressure, temperature, ok)
  if (.not. ok) error stop "Failed to evaluate right HLLC density"

  ! Equal states must reduce exactly to the physical flux.
  velocity = 75.0_dp
  call set_primitive(rho_l, velocity, pressure, yl, ql)
  call reactive_primitive_to_conserved(species, ql, ul, tl, cl, ok)
  if (.not. ok) error stop "Failed to build equal HLLC state"
  call reactive_hllc_flux_x(species, ul, ul, tl, tl, flux, ok)
  if (.not. ok) error stop "Equal-state HLLC flux failed"
  call reactive_physical_flux_x( &
    species, ul, tl, reference, dummy_t, dummy_c, primitive, ok)
  if (.not. ok) error stop "Equal-state physical flux failed"
  call assert_vector_close(flux, reference, 2.0e-12_dp, &
    "equal-state physical flux")
  call reactive_pelec_flux_x(species, ul, ul, tl, tl, pelec, ok)
  if (.not. ok) error stop "Equal-state reactive PeleC flux failed"
  call assert_vector_close(pelec, reference, 2.0e-12_dp, &
    "equal-state PeleC physical flux")

  ! A stationary material contact has zero mass/energy/species flux and p momentum.
  velocity = 0.0_dp
  call set_primitive(rho_l, velocity, pressure, yl, ql)
  call set_primitive(rho_r, velocity, pressure, yr, qr)
  call reactive_primitive_to_conserved(species, ql, ul, tl, cl, ok)
  if (.not. ok) error stop "Failed to build stationary left contact"
  call reactive_primitive_to_conserved(species, qr, ur, tr, cr, ok)
  if (.not. ok) error stop "Failed to build stationary right contact"
  call reactive_hllc_flux_x(species, ul, ur, tl, tr, flux, ok)
  if (.not. ok) error stop "Stationary-contact HLLC flux failed"
  call assert_close(flux(irho), 0.0_dp, 2.0e-12_dp, "stationary mass flux")
  call assert_close(flux(imx), pressure, 2.0e-12_dp, &
    "stationary momentum flux")
  call assert_close(flux(imy), 0.0_dp, 2.0e-12_dp, "stationary y flux")
  call assert_close(flux(imz), 0.0_dp, 2.0e-12_dp, "stationary z flux")
  call assert_close(flux(iet), 0.0_dp, 2.0e-12_dp, "stationary energy flux")
  call assert_close(sum(abs(flux(6:12))), 0.0_dp, 2.0e-12_dp, &
    "stationary species flux")
  call reactive_pelec_flux_x(species, ul, ur, tl, tr, pelec, ok)
  if (.not. ok) error stop "Stationary-contact reactive PeleC flux failed"
  call assert_close(pelec(irho), 0.0_dp, 2.0e-12_dp, &
    "PeleC stationary mass flux")
  call assert_close(pelec(imx), pressure, 2.0e-12_dp, &
    "PeleC stationary momentum flux")
  call assert_close(pelec(iet), 0.0_dp, 2.0e-12_dp, &
    "PeleC stationary energy flux")
  call assert_close(sum(abs(pelec(6:12))), 0.0_dp, 2.0e-12_dp, &
    "PeleC stationary species flux")
  call reactive_rusanov_flux_x(species, ul, ur, tl, tr, rusanov, ok)
  if (.not. ok) error stop "Stationary-contact Rusanov flux failed"
  if (abs(rusanov(irho)) <= 1.0e-3_dp) then
    error stop "Rusanov baseline unexpectedly resolved stationary contact"
  end if

  ! A moving material contact is exactly upwinded by HLLC.
  velocity = 100.0_dp
  call set_primitive(rho_l, velocity, pressure, yl, ql)
  call set_primitive(rho_r, velocity, pressure, yr, qr)
  call reactive_primitive_to_conserved(species, ql, ul, tl, cl, ok)
  if (.not. ok) error stop "Failed to build moving left contact"
  call reactive_primitive_to_conserved(species, qr, ur, tr, cr, ok)
  if (.not. ok) error stop "Failed to build moving right contact"
  call reactive_hllc_flux_x(species, ul, ur, tl, tr, flux, ok)
  if (.not. ok) error stop "Moving-contact HLLC flux failed"
  call reactive_physical_flux_x( &
    species, ul, tl, reference, dummy_t, dummy_c, primitive, ok)
  if (.not. ok) error stop "Moving-contact physical flux failed"
  call assert_vector_close(flux, reference, 5.0e-12_dp, &
    "moving-contact upwind flux")
  call reactive_pelec_flux_x(species, ul, ur, tl, tr, pelec, ok)
  if (.not. ok) error stop "Moving-contact reactive PeleC flux failed"
  call assert_vector_close(pelec, reference, 5.0e-12_dp, &
    "moving-contact PeleC upwind flux")
  call assert_close(sum(flux(6:12)), flux(irho), 2.0e-12_dp, &
    "HLLC species closure")

  ! A finite pressure jump exercises the acoustic-star paths away from contacts.
  xl = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]
  call mass_fractions_from_mole_fractions(species, xl, yl, ok)
  if (.not. ok) error stop "Failed to build acoustic composition"
  rho_l = mixture_density(species, yl, 2.0_dp * pressure, 1200.0_dp, ok)
  if (.not. ok) error stop "Failed to evaluate acoustic left density"
  rho_r = mixture_density(species, yl, pressure, 900.0_dp, ok)
  if (.not. ok) error stop "Failed to evaluate acoustic right density"
  call set_primitive(rho_l, 0.0_dp, 2.0_dp * pressure, yl, ql)
  call set_primitive(rho_r, 0.0_dp, pressure, yl, qr)
  call reactive_primitive_to_conserved(species, ql, ul, tl, cl, ok)
  if (.not. ok) error stop "Failed to build acoustic left state"
  call reactive_primitive_to_conserved(species, qr, ur, tr, cr, ok)
  if (.not. ok) error stop "Failed to build acoustic right state"
  call reactive_hllc_flux_x(species, ul, ur, tl, tr, flux, ok)
  if (.not. ok) error stop "Acoustic HLLC flux failed"
  if (flux(irho) <= 0.0_dp .or. flux(imx) <= pressure .or. &
      flux(iet) <= 0.0_dp) then
    error stop "Acoustic HLLC flux has the wrong propagation direction"
  end if
  call assert_close(sum(flux(6:12)), flux(irho), 2.0e-12_dp, &
    "acoustic HLLC species closure")
  call reactive_pelec_flux_x(species, ul, ur, tl, tr, pelec, ok)
  if (.not. ok) error stop "Acoustic reactive PeleC flux failed"
  if (pelec(irho) <= 0.0_dp .or. pelec(imx) <= pressure .or. &
      pelec(iet) <= 0.0_dp) then
    error stop "Acoustic PeleC flux has the wrong propagation direction"
  end if
  call assert_close(sum(pelec(6:12)), pelec(irho), 2.0e-12_dp, &
    "acoustic PeleC species closure")

  call compute_reactive_riemann_flux_x( &
    species, ul, ur, tl, tr, "hllc", flux, ok)
  if (.not. ok) error stop "HLLC dispatch failed"
  call compute_reactive_riemann_flux_x( &
    species, ul, ur, tl, tr, "pelec", flux, ok)
  if (.not. ok) error stop "Reactive PeleC dispatch failed"
  call compute_reactive_riemann_flux_x( &
    species, ul, ur, tl, tr, "not-a-solver", flux, ok)
  if (ok) error stop "Unknown reactive Riemann solver was accepted"

  error = abs(sum(yl) - 1.0_dp) + abs(sum(yr) - 1.0_dp)
  if (error > 2.0e-14_dp) error stop "HLLC test composition closure failed"
  write(*, '(a)') "test_reactive_riemann: PASS"

contains

  subroutine set_primitive(rho, u, p, y, q)
    real(dp), intent(in) :: rho, u, p, y(:)
    real(dp), intent(out) :: q(:)
    integer :: species_index

    q = 0.0_dp
    q(1:5) = [rho, u, 0.0_dp, 0.0_dp, p]
    do species_index = 1, size(y)
      q(reactive_mass_fraction_component(species_index)) = y(species_index)
    end do
  end subroutine set_primitive

  subroutine assert_close(actual, expected, relative_tolerance, label)
    real(dp), intent(in) :: actual, expected, relative_tolerance
    character(len=*), intent(in) :: label
    real(dp) :: scaled_error

    scaled_error = abs(actual - expected) / max(1.0_dp, abs(expected))
    if (scaled_error > relative_tolerance) then
      write(*, '(a,3(1x,es24.16))') &
        trim(label), actual, expected, scaled_error
      error stop "Reactive Riemann scalar mismatch"
    end if
  end subroutine assert_close

  subroutine assert_vector_close(actual, expected, relative_tolerance, label)
    real(dp), intent(in) :: actual(:), expected(:), relative_tolerance
    character(len=*), intent(in) :: label
    real(dp) :: scaled_error

    scaled_error = maxval(abs(actual - expected) / max(1.0_dp, abs(expected)))
    if (scaled_error > relative_tolerance) then
      write(*, '(a,1x,es24.16)') trim(label), scaled_error
      error stop "Reactive Riemann vector mismatch"
    end if
  end subroutine assert_vector_close

end program test_reactive_riemann
