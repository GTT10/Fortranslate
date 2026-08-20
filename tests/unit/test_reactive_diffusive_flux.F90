program test_reactive_diffusive_flux
  use precision_mod, only: dp
  use state_indices_mod, only: irho, imx, imy, imz, iet
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: &
    mass_fractions_from_mole_fractions, mixture_density
  use transport_database_mod, only: &
    gas_transport_species, load_h2o2_elementary_transport
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_species_component, reactive_primitive_to_conserved, &
    reactive_diffusive_flux_x
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  real(dp), allocatable :: qleft(:), qright(:), uleft(:), uright(:), flux(:)
  real(dp) :: xbase(7), xleft(7), xright(7), yleft(7), yright(7)
  real(dp) :: tleft, tright, sound_speed, mu, lambda, diffusion(7), dx
  logical :: ok
  integer :: k, nvar, nprim

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "thermodynamics load")
  call load_h2o2_elementary_transport(transport, ok)
  call require(ok, "transport database loads")
  nvar = reactive_nvar(size(species))
  nprim = reactive_nprim(size(species))
  allocate(qleft(nprim), qright(nprim), uleft(nvar), uright(nvar), flux(nvar))
  dx = 2.0e-4_dp
  xbase = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, 1.0e-5_dp, &
    0.0_dp, 0.55643_dp]
  xbase = xbase / sum(xbase)
  call mass_fractions_from_mole_fractions(species, xbase, yleft, ok)
  call require(ok, "base composition converts")
  yright = yleft

  call make_state(1000.0_dp, 101325.0_dp, 0.0_dp, 0.0_dp, yleft, &
    qleft, uleft, tleft)
  qright = qleft
  uright = uleft
  tright = tleft
  call reactive_diffusive_flux_x( &
    species, transport, uleft, uright, tleft, tright, dx, .true., .true., &
    .true., .true., flux, ok, mu, lambda, diffusion)
  call require(ok, "equal-state diffusive flux evaluates")
  call require(maxval(abs(flux)) < 1.0e-18_dp, &
    "equal state has zero diffusive flux")

  call make_state(1000.0_dp, 101325.0_dp, 0.0_dp, 0.0_dp, yleft, &
    qleft, uleft, tleft)
  call make_state(1000.0_dp, 101325.0_dp, 0.0_dp, 1.0_dp, yright, &
    qright, uright, tright)
  call reactive_diffusive_flux_x( &
    species, transport, uleft, uright, tleft, tright, dx, .true., .false., &
    .false., .false., flux, ok, mu, lambda, diffusion)
  call require(ok, "viscous face flux evaluates")
  call require_close(flux(imy), -mu / dx, 2.0e-12_dp, &
    "transverse viscous momentum flux")
  call require_close(flux(iet), -0.5_dp * mu / dx, 2.0e-12_dp, &
    "viscous work is included in energy flux")
  call require(abs(flux(irho)) < 1.0e-30_dp, &
    "viscosity does not diffuse total mass")

  call make_state(900.0_dp, 101325.0_dp, 0.0_dp, 0.0_dp, yleft, &
    qleft, uleft, tleft)
  call make_state(1100.0_dp, 101325.0_dp, 0.0_dp, 0.0_dp, yright, &
    qright, uright, tright)
  call reactive_diffusive_flux_x( &
    species, transport, uleft, uright, tleft, tright, dx, .false., .true., &
    .false., .false., flux, ok, mu, lambda, diffusion)
  call require(ok, "conductive face flux evaluates")
  call require_close(flux(iet), -lambda * (tright - tleft) / dx, &
    2.0e-12_dp, "Fourier heat flux")
  call require(flux(iet) < 0.0_dp, &
    "heat flux points from the hotter right cell to the left")

  xleft = xbase
  xright = xbase
  xleft(1) = xleft(1) - 0.03_dp
  xleft(7) = xleft(7) + 0.03_dp
  xright(1) = xright(1) + 0.03_dp
  xright(7) = xright(7) - 0.03_dp
  call mass_fractions_from_mole_fractions(species, xleft, yleft, ok)
  call require(ok, "left diffusion composition converts")
  call mass_fractions_from_mole_fractions(species, xright, yright, ok)
  call require(ok, "right diffusion composition converts")
  call make_state(1000.0_dp, 101325.0_dp, 0.0_dp, 0.0_dp, yleft, &
    qleft, uleft, tleft)
  call make_state(1000.0_dp, 101325.0_dp, 0.0_dp, 0.0_dp, yright, &
    qright, uright, tright)
  call reactive_diffusive_flux_x( &
    species, transport, uleft, uright, tleft, tright, dx, .false., .false., &
    .true., .true., flux, ok, mu, lambda, diffusion)
  call require(ok, "species-diffusion face flux evaluates")
  call require(abs(sum(flux(reactive_species_component(1): &
    reactive_species_component(7)))) < 2.0e-15_dp, &
    "correction velocity closes species mass flux")
  call require(flux(reactive_species_component(1)) < 0.0_dp, &
    "H2 diffuses down its mole-fraction gradient")
  call require(abs(flux(irho)) < 1.0e-30_dp, &
    "species diffusion does not create total-mass flux")
  call require(abs(flux(iet)) > 1.0e-8_dp, &
    "species enthalpy transport contributes to energy flux")

  call make_state(1000.0_dp, 0.98_dp * 101325.0_dp, 0.0_dp, 0.0_dp, &
    yleft, qleft, uleft, tleft)
  call make_state(1000.0_dp, 1.02_dp * 101325.0_dp, 0.0_dp, 0.0_dp, &
    yleft, qright, uright, tright)
  call reactive_diffusive_flux_x( &
    species, transport, uleft, uright, tleft, tright, dx, .false., .false., &
    .true., .false., flux, ok, mu, lambda, diffusion)
  call require(ok, "pressure-gradient reference flux evaluates")
  call require(maxval(abs(flux(reactive_species_component(1): &
    reactive_species_component(7)))) < 1.0e-18_dp, &
    "species flux is zero without barodiffusion")
  call reactive_diffusive_flux_x( &
    species, transport, uleft, uright, tleft, tright, dx, .false., .false., &
    .true., .true., flux, ok, mu, lambda, diffusion)
  call require(ok, "barodiffusive face flux evaluates")
  call require(maxval(abs(flux(reactive_species_component(1): &
    reactive_species_component(7)))) > 1.0e-8_dp, &
    "pressure gradient activates barodiffusion")
  call require(abs(sum(flux(reactive_species_component(1): &
    reactive_species_component(7)))) < 2.0e-15_dp, &
    "barodiffusive correction velocity closes mass flux")

contains

  subroutine make_state(temperature, pressure, u, v, y, q, conserved, &
      recovered_temperature)
    real(dp), intent(in) :: temperature, pressure, u, v, y(:)
    real(dp), intent(out) :: q(:), conserved(:), recovered_temperature
    real(dp) :: rho, c
    logical :: local_ok
    integer :: n

    rho = mixture_density(species, y, pressure, temperature, local_ok)
    call require(local_ok, "density evaluates")
    q(1:5) = [rho, u, v, 0.0_dp, pressure]
    do n = 1, size(species)
      q(reactive_mass_fraction_component(n)) = y(n)
    end do
    call reactive_primitive_to_conserved( &
      species, q, conserved, recovered_temperature, c, local_ok)
    call require(local_ok, "primitive state converts")
  end subroutine make_state

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') "FAILED: " // trim(message)
      error stop 1
    end if
  end subroutine require

  subroutine require_close(actual, expected, relative_tolerance, message)
    real(dp), intent(in) :: actual, expected, relative_tolerance
    character(len=*), intent(in) :: message
    real(dp) :: error
    error = abs(actual - expected) / max(1.0e-30_dp, abs(expected))
    call require(error <= relative_tolerance, message)
  end subroutine require_close

end program test_reactive_diffusive_flux
