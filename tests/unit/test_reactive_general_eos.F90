program test_reactive_general_eos
  use precision_mod, only: dp
  use state_indices_mod, only: irho, imx, imy, imz, iet
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_primitive_to_conserved, reactive_conserved_to_primitive, &
    reactive_rusanov_flux_x, reactive_difference_to_characteristics, &
    reactive_characteristics_to_difference, trace_reactive_characteristics
  implicit none

  type(nasa7_species), allocatable :: species(:)
  real(dp), allocatable :: q(:), recovered(:), conserved(:), flux(:)
  real(dp), allocatable :: difference(:), reconstructed(:), left(:), right(:)
  real(dp) :: mole_fractions(7), mass_fractions(7)
  real(dp) :: temperature, recovered_temperature, sound_speed, recovered_c
  real(dp) :: characteristic(5), tolerance
  logical :: ok
  integer :: k

  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok) error stop "Failed to load reactive thermodynamics"
  allocate(q(reactive_nprim(7)), recovered(reactive_nprim(7)))
  allocate(conserved(reactive_nvar(7)), flux(reactive_nvar(7)))
  allocate(difference(reactive_nprim(7)), reconstructed(reactive_nprim(7)))
  allocate(left(reactive_nprim(7)), right(reactive_nprim(7)))

  mole_fractions = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
    1.0e-5_dp, 0.0_dp, 0.55643_dp]
  call mass_fractions_from_mole_fractions( &
    species, mole_fractions, mass_fractions, ok)
  if (.not. ok) error stop "Failed to convert reactive composition"

  q(1:5) = [0.45_dp, 120.0_dp, -4.0_dp, 2.5_dp, 2.0e5_dp]
  do k = 1, 7
    q(reactive_mass_fraction_component(k)) = mass_fractions(k)
  end do
  call reactive_primitive_to_conserved( &
    species, q, conserved, temperature, sound_speed, ok)
  if (.not. ok) error stop "Primitive-to-conserved conversion failed"
  call reactive_conserved_to_primitive( &
    species, conserved, temperature, recovered, recovered_temperature, &
    recovered_c, ok)
  if (.not. ok) error stop "Conserved-to-primitive conversion failed"
  call assert_vector_close(recovered, q, 2.0e-11_dp, "state round trip")
  call assert_close(recovered_temperature, temperature, 2.0e-11_dp, &
    "temperature round trip")
  call assert_close(recovered_c, sound_speed, 2.0e-11_dp, &
    "sound speed round trip")

  call reactive_rusanov_flux_x( &
    species, conserved, conserved, temperature, temperature, flux, ok)
  if (.not. ok) error stop "Equal-state reactive flux failed"
  call assert_close(flux(irho), q(1) * q(2), 2.0e-12_dp, "mass flux")
  call assert_close(flux(imx), q(1) * q(2)**2 + q(5), &
    2.0e-12_dp, "x-momentum flux")
  call assert_close(flux(imy), q(1) * q(2) * q(3), &
    2.0e-12_dp, "y-momentum flux")
  call assert_close(flux(imz), q(1) * q(2) * q(4), &
    2.0e-12_dp, "z-momentum flux")
  call assert_close(flux(iet), (conserved(iet) + q(5)) * q(2), &
    2.0e-12_dp, "energy flux")
  call assert_close(sum(flux(6:12)), flux(irho), 2.0e-12_dp, &
    "species-flux closure")

  difference = 0.0_dp
  difference(1:5) = [0.02_dp, -3.0_dp, 0.4_dp, -0.2_dp, 2500.0_dp]
  call reactive_difference_to_characteristics( &
    q, difference, sound_speed, characteristic, ok)
  if (.not. ok) error stop "Characteristic projection failed"
  call reactive_characteristics_to_difference( &
    q, characteristic, sound_speed, reconstructed, ok)
  if (.not. ok) error stop "Characteristic inverse failed"
  call assert_vector_close(reconstructed(1:5), difference(1:5), &
    2.0e-12_dp, "characteristic round trip")

  do k = 1, 7
    difference(reactive_mass_fraction_component(k)) = &
      1.0e-4_dp * real(k - 4, dp)
  end do
  call trace_reactive_characteristics( &
    q, difference, sound_speed, 0.0_dp, left, right, ok)
  if (.not. ok) error stop "Zero-Courant characteristic trace failed"
  call assert_vector_close(left, q - 0.5_dp * difference, &
    2.0e-12_dp, "left zero-Courant state")
  call assert_vector_close(right, q + 0.5_dp * difference, &
    2.0e-12_dp, "right zero-Courant state")

  tolerance = abs(sum(conserved(6:12)) - conserved(irho))
  if (tolerance > 2.0e-13_dp) error stop "Conserved species closure failed"
  write(*, '(a)') "test_reactive_general_eos: PASS"

contains

  subroutine assert_close(actual, expected, relative_tolerance, label)
    real(dp), intent(in) :: actual, expected, relative_tolerance
    character(len=*), intent(in) :: label
    real(dp) :: error
    error = abs(actual - expected) / max(1.0_dp, abs(expected))
    if (error > relative_tolerance) then
      write(*, '(a,3(1x,es24.16))') trim(label), actual, expected, error
      error stop "Reactive general-EOS mismatch"
    end if
  end subroutine assert_close

  subroutine assert_vector_close(actual, expected, relative_tolerance, label)
    real(dp), intent(in) :: actual(:), expected(:), relative_tolerance
    character(len=*), intent(in) :: label
    real(dp) :: error
    error = maxval(abs(actual - expected) / max(1.0_dp, abs(expected)))
    if (error > relative_tolerance) then
      write(*, '(a,1x,es24.16)') trim(label), error
      error stop "Reactive vector mismatch"
    end if
  end subroutine assert_vector_close

end program test_reactive_general_eos
