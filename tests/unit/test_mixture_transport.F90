program test_mixture_transport
  use, intrinsic :: ieee_arithmetic, only: &
    ieee_is_finite, ieee_positive_inf, ieee_quiet_nan, ieee_value
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use gas_transport_mod, only: &
    gas_transport_species, valid_gas_transport_species, &
    compatible_transport_database
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions
  use transport_database_mod, only: load_h2o2_elementary_transport
  use mixture_transport_mod, only: &
    collision_integral_diffusion, collision_integral_viscosity, &
    pure_species_viscosities, pure_species_thermal_conductivities, &
    mixture_viscosity_wilke, mixture_thermal_conductivity_mathur, &
    binary_diffusion_coefficients, mixture_averaged_diffusion_coefficients, &
    mixture_transport_coefficients, standard_atmosphere
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  real(dp) :: x(7), y(7), mu, lambda, diffusion(7)
  real(dp) :: pure_mu(7), binary_1(7, 7), binary_2(7, 7)
  logical :: ok

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "thermodynamic database loads")
  call load_h2o2_elementary_transport(transport, ok)
  call require(ok, "transport database loads")
  call require(compatible_transport_database(species, transport), &
    "transport and thermodynamic species ordering agrees")

  call pure_species_viscosities(species, transport, 1000.0_dp, pure_mu, ok)
  call require(ok, "pure-species viscosities evaluate")
  call require_close(pure_mu(1), 1.9684274672351392e-5_dp, 2.0e-12_dp, &
    "H2 viscosity at 1000 K")
  call require_close(pure_mu(4), 4.7888632094172495e-5_dp, 2.0e-12_dp, &
    "O2 viscosity at 1000 K")
  call require_close(pure_mu(7), 4.146512036043256e-5_dp, 2.0e-12_dp, &
    "N2 viscosity at 1000 K")

  x = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, 1.0e-5_dp, &
    0.0_dp, 0.55643_dp]
  x = x / sum(x)
  call mass_fractions_from_mole_fractions(species, x, y, ok)
  call require(ok, "reference composition converts")
  call mixture_transport_coefficients( &
    species, transport, y, 1000.0_dp, standard_atmosphere, mu, lambda, &
    diffusion, ok)
  call require(ok, "mixture transport evaluates")
  call require_close(mu, 4.1983389803134420e-5_dp, 2.0e-12_dp, &
    "Wilke mixture viscosity")
  call require_close(lambda, 1.2428859732637745e-1_dp, 2.0e-12_dp, &
    "Mathur mixture conductivity")
  call require_close(diffusion(1), 8.1296448205293609e-4_dp, 2.0e-12_dp, &
    "H2 mixture-averaged diffusion")
  call require_close(diffusion(4), 1.9833619301130086e-4_dp, 2.0e-12_dp, &
    "O2 mixture-averaged diffusion")
  call require_close(diffusion(7), 1.8023444527551745e-4_dp, 2.0e-12_dp, &
    "N2 mixture-averaged diffusion")

  call binary_diffusion_coefficients( &
    species, transport, 1000.0_dp, standard_atmosphere, binary_1, ok)
  call require(ok, "binary diffusion evaluates")
  call binary_diffusion_coefficients( &
    species, transport, 1000.0_dp, 2.0_dp * standard_atmosphere, binary_2, ok)
  call require(ok, "binary diffusion evaluates at doubled pressure")
  call require(maxval(abs(binary_1 - transpose(binary_1))) < 1.0e-18_dp, &
    "binary diffusion matrix is symmetric")
  call require(maxval(abs(diagonal(binary_1))) < 1.0e-30_dp, &
    "binary diffusion diagonal is zero")
  call require_close(binary_2(1, 4), 0.5_dp * binary_1(1, 4), &
    2.0e-13_dp, "binary diffusion scales inversely with pressure")

  call check_invalid_inputs(species, transport, y, ok)
  call require(ok, "invalid transport inputs fail closed")
  write(*, '(a)') "test_mixture_transport: PASS"

contains

  subroutine check_invalid_inputs(species, transport, mass_fractions, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: mass_fractions(:)
    logical, intent(out) :: ok

    type(gas_transport_species) :: bad_transport(size(transport))
    real(dp) :: bad_mass_fractions(size(mass_fractions))
    real(dp) :: nan_value, positive_inf
    real(dp) :: pure_values(size(species)), binary_values(size(species), size(species))
    real(dp) :: diffusion_values(size(species))
    real(dp) :: viscosity, conductivity
    logical :: local_ok

    ok = .false.
    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    positive_inf = ieee_value(0.0_dp, ieee_positive_inf)
    call require(.not. ieee_is_finite(nan_value), "quiet NaN is available")
    call require(.not. ieee_is_finite(positive_inf), "positive infinity is available")

    bad_transport = transport
    bad_transport(1)%well_depth = nan_value
    call require(.not. valid_gas_transport_species(bad_transport(1)), &
      "NaN transport record is rejected")
    bad_transport = transport
    bad_transport(1)%diameter = positive_inf
    call require(.not. valid_gas_transport_species(bad_transport(1)), &
      "+Inf transport record is rejected")

    call require(collision_integral_viscosity(nan_value) < 0.0_dp, &
      "NaN viscosity collision integral is rejected")
    call require(collision_integral_viscosity(positive_inf) < 0.0_dp, &
      "+Inf viscosity collision integral is rejected")
    call require(collision_integral_diffusion(nan_value) < 0.0_dp, &
      "NaN diffusion collision integral is rejected")
    call require(collision_integral_diffusion(positive_inf) < 0.0_dp, &
      "+Inf diffusion collision integral is rejected")
    viscosity = collision_integral_viscosity(huge(1.0_dp))
    call require(ieee_is_finite(viscosity) .and. viscosity > 0.0_dp, &
      "huge finite viscosity collision argument remains finite")
    viscosity = collision_integral_diffusion(huge(1.0_dp))
    call require(ieee_is_finite(viscosity) .and. viscosity > 0.0_dp, &
      "huge finite diffusion collision argument remains finite")

    pure_values = 1.0_dp
    call pure_species_viscosities( &
      species, transport, nan_value, pure_values, local_ok)
    call require(.not. local_ok .and. all(pure_values == 0.0_dp), &
      "NaN temperature viscosity output is zero")
    pure_values = 1.0_dp
    call pure_species_viscosities( &
      species, transport, positive_inf, pure_values, local_ok)
    call require(.not. local_ok .and. all(pure_values == 0.0_dp), &
      "+Inf temperature viscosity output is zero")

    pure_values = 1.0_dp
    call pure_species_thermal_conductivities( &
      species, transport, nan_value, pure_values, local_ok)
    call require(.not. local_ok .and. all(pure_values == 0.0_dp), &
      "NaN temperature conductivity output is zero")
    pure_values = 1.0_dp
    call pure_species_thermal_conductivities( &
      species, transport, positive_inf, pure_values, local_ok)
    call require(.not. local_ok .and. all(pure_values == 0.0_dp), &
      "+Inf temperature conductivity output is zero")

    viscosity = 1.0_dp
    call mixture_viscosity_wilke( &
      species, transport, mass_fractions, nan_value, viscosity, local_ok)
    call require(.not. local_ok .and. viscosity == 0.0_dp, &
      "NaN mixture temperature viscosity output is zero")
    viscosity = 1.0_dp
    call mixture_viscosity_wilke( &
      species, transport, mass_fractions, positive_inf, viscosity, local_ok)
    call require(.not. local_ok .and. viscosity == 0.0_dp, &
      "+Inf mixture temperature viscosity output is zero")

    conductivity = 1.0_dp
    call mixture_thermal_conductivity_mathur( &
      species, transport, mass_fractions, nan_value, conductivity, local_ok)
    call require(.not. local_ok .and. conductivity == 0.0_dp, &
      "NaN mixture temperature conductivity output is zero")
    conductivity = 1.0_dp
    call mixture_thermal_conductivity_mathur( &
      species, transport, mass_fractions, positive_inf, conductivity, local_ok)
    call require(.not. local_ok .and. conductivity == 0.0_dp, &
      "+Inf mixture temperature conductivity output is zero")

    binary_values = 1.0_dp
    call binary_diffusion_coefficients( &
      species, transport, nan_value, standard_atmosphere, binary_values, local_ok)
    call require(.not. local_ok .and. all(binary_values == 0.0_dp), &
      "NaN temperature binary output is zero")
    binary_values = 1.0_dp
    call binary_diffusion_coefficients( &
      species, transport, 1000.0_dp, positive_inf, binary_values, local_ok)
    call require(.not. local_ok .and. all(binary_values == 0.0_dp), &
      "+Inf pressure binary output is zero")
    binary_values = 1.0_dp
    call binary_diffusion_coefficients( &
      species, transport, 1000.0_dp, nan_value, binary_values, local_ok)
    call require(.not. local_ok .and. all(binary_values == 0.0_dp), &
      "NaN pressure binary output is zero")

    diffusion_values = 1.0_dp
    call mixture_averaged_diffusion_coefficients( &
      species, transport, mass_fractions, 1000.0_dp, positive_inf, &
      diffusion_values, local_ok)
    call require(.not. local_ok .and. all(diffusion_values == 0.0_dp), &
      "+Inf pressure mixture diffusion output is zero")
    diffusion_values = 1.0_dp
    call mixture_averaged_diffusion_coefficients( &
      species, transport, mass_fractions, 1000.0_dp, nan_value, &
      diffusion_values, local_ok)
    call require(.not. local_ok .and. all(diffusion_values == 0.0_dp), &
      "NaN pressure mixture diffusion output is zero")

    bad_mass_fractions = mass_fractions
    bad_mass_fractions(1) = nan_value
    viscosity = 1.0_dp
    call mixture_viscosity_wilke( &
      species, transport, bad_mass_fractions, 1000.0_dp, viscosity, local_ok)
    call require(.not. local_ok .and. viscosity == 0.0_dp, &
      "NaN composition viscosity output is zero")
    bad_mass_fractions = mass_fractions
    bad_mass_fractions(1) = positive_inf
    conductivity = 1.0_dp
    call mixture_thermal_conductivity_mathur( &
      species, transport, bad_mass_fractions, 1000.0_dp, conductivity, local_ok)
    call require(.not. local_ok .and. conductivity == 0.0_dp, &
      "+Inf composition conductivity output is zero")
    diffusion_values = 1.0_dp
    call mixture_averaged_diffusion_coefficients( &
      species, transport, bad_mass_fractions, 1000.0_dp, standard_atmosphere, &
      diffusion_values, local_ok)
    call require(.not. local_ok .and. all(diffusion_values == 0.0_dp), &
      "+Inf composition diffusion output is zero")

    bad_transport = transport
    bad_transport(1)%diameter = huge(1.0_dp)
    pure_values = 1.0_dp
    call pure_species_viscosities( &
      species, bad_transport, 1000.0_dp, pure_values, local_ok)
    call require(.not. local_ok .and. all(pure_values == 0.0_dp), &
      "overflowing transport intermediate is rejected")
    binary_values = 1.0_dp
    call binary_diffusion_coefficients( &
      species, bad_transport, 1000.0_dp, standard_atmosphere, binary_values, local_ok)
    call require(.not. local_ok .and. all(binary_values == 0.0_dp), &
      "overflowing binary intermediate is rejected")

    bad_transport = transport
    bad_transport(1)%well_depth = tiny(1.0_dp)
    pure_values = 1.0_dp
    call pure_species_viscosities( &
      species, bad_transport, 1000.0_dp, pure_values, local_ok)
    call require(.not. local_ok .and. all(pure_values == 0.0_dp), &
      "overflowing reduced temperature is rejected")

    viscosity = 1.0_dp
    conductivity = 1.0_dp
    diffusion_values = 1.0_dp
    call mixture_transport_coefficients( &
      species, transport, mass_fractions, nan_value, standard_atmosphere, &
      viscosity, conductivity, diffusion_values, local_ok)
    call require(.not. local_ok .and. viscosity == 0.0_dp .and. &
      conductivity == 0.0_dp .and. all(diffusion_values == 0.0_dp), &
      "invalid combined transport outputs are zero")

    ok = .true.
  end subroutine check_invalid_inputs

  pure function diagonal(matrix) result(values)
    real(dp), intent(in) :: matrix(:, :)
    real(dp) :: values(min(size(matrix, 1), size(matrix, 2)))
    integer :: i
    do i = 1, size(values)
      values(i) = matrix(i, i)
    end do
  end function diagonal

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

end program test_mixture_transport
