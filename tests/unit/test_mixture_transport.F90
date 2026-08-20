program test_mixture_transport
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions
  use transport_database_mod, only: &
    gas_transport_species, load_h2o2_elementary_transport, &
    compatible_transport_database
  use mixture_transport_mod, only: &
    pure_species_viscosities, binary_diffusion_coefficients, &
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

contains

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
