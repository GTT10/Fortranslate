program test_reactive_general_eos
  use, intrinsic :: ieee_arithmetic, only: &
    ieee_is_finite, ieee_positive_inf, ieee_quiet_nan, ieee_value
  use precision_mod, only: dp
  use state_indices_mod, only: irho, imx, imy, imz, iet
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_species_component, &
    reactive_mass_fraction_component, &
    reactive_primitive_to_conserved, reactive_conserved_to_primitive, &
    reactive_rusanov_flux_x, reactive_difference_to_characteristics, &
    reactive_characteristics_to_difference, trace_reactive_characteristics
  use reactive_2d_mod, only: compute_reactive_cfl_timestep_2d
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
  call check_transform_failure_contract(species, q, conserved, temperature)
  call check_2d_timestep_failure_contract(species, conserved, temperature)
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

  subroutine require(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label

    if (.not. condition) error stop trim(label)
  end subroutine require

  subroutine check_transform_failure_contract( &
      active_species, valid_primitive, valid_conserved, temperature_guess)
    type(nasa7_species), intent(in) :: active_species(:)
    real(dp), intent(in) :: valid_primitive(:), valid_conserved(:)
    real(dp), intent(in) :: temperature_guess

    real(dp), allocatable :: bad_primitive(:), bad_conserved(:)
    real(dp) :: nan_value, positive_inf

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    positive_inf = ieee_value(0.0_dp, ieee_positive_inf)
    call require(.not. ieee_is_finite(nan_value), "quiet NaN unavailable")
    call require(.not. ieee_is_finite(positive_inf), &
      "positive infinity unavailable")
    allocate(bad_primitive(size(valid_primitive)))
    allocate(bad_conserved(size(valid_conserved)))

    bad_primitive = valid_primitive
    bad_primitive(1) = nan_value
    call require_primitive_rejected( &
      active_species, bad_primitive, "primitive NaN density")
    bad_primitive = valid_primitive
    bad_primitive(5) = positive_inf
    call require_primitive_rejected( &
      active_species, bad_primitive, "primitive infinite pressure")
    bad_primitive = valid_primitive
    bad_primitive(2) = huge(1.0_dp)
    call require_primitive_rejected( &
      active_species, bad_primitive, "primitive extreme velocity")
    bad_primitive = valid_primitive
    bad_primitive(reactive_mass_fraction_component(1)) = nan_value
    call require_primitive_rejected( &
      active_species, bad_primitive, "primitive NaN composition")
    bad_primitive = valid_primitive
    bad_primitive(reactive_mass_fraction_component(1)) = positive_inf
    call require_primitive_rejected( &
      active_species, bad_primitive, "primitive infinite composition")

    bad_conserved = valid_conserved
    bad_conserved(irho) = nan_value
    call require_conserved_rejected( &
      active_species, bad_conserved, temperature_guess, &
      "conserved NaN density")
    bad_conserved = valid_conserved
    bad_conserved(iet) = positive_inf
    call require_conserved_rejected( &
      active_species, bad_conserved, temperature_guess, &
      "conserved infinite energy")
    bad_conserved = valid_conserved
    bad_conserved(imx) = huge(1.0_dp)
    call require_conserved_rejected( &
      active_species, bad_conserved, temperature_guess, &
      "conserved extreme momentum")
    bad_conserved = valid_conserved
    bad_conserved(reactive_species_component(1)) = nan_value
    call require_conserved_rejected( &
      active_species, bad_conserved, temperature_guess, &
      "conserved NaN species")
    bad_conserved = valid_conserved
    bad_conserved(reactive_species_component(1)) = positive_inf
    call require_conserved_rejected( &
      active_species, bad_conserved, temperature_guess, &
      "conserved infinite species")
    bad_conserved = valid_conserved
    bad_conserved(reactive_species_component(2)) = huge(1.0_dp)
    call require_conserved_rejected( &
      active_species, bad_conserved, temperature_guess, &
      "conserved extreme species density")
    bad_conserved = valid_conserved
    call require_conserved_rejected( &
      active_species, bad_conserved, nan_value, &
      "conserved NaN temperature guess")
    call require_conserved_rejected( &
      active_species, bad_conserved, positive_inf, &
      "conserved infinite temperature guess")
  end subroutine check_transform_failure_contract

  subroutine require_primitive_rejected(active_species, input, label)
    type(nasa7_species), intent(in) :: active_species(:)
    real(dp), intent(in) :: input(:)
    character(len=*), intent(in) :: label

    real(dp), allocatable :: output(:)
    real(dp) :: output_temperature, output_sound_speed
    logical :: local_ok

    allocate(output(reactive_nvar(size(active_species))))
    output = 7.0_dp
    output_temperature = 7.0_dp
    output_sound_speed = 7.0_dp
    call reactive_primitive_to_conserved( &
      active_species, input, output, output_temperature, output_sound_speed, &
      local_ok)
    call require(.not. local_ok, trim(label) // " accepted")
    call require(all(output == 0.0_dp), trim(label) // " output changed")
    call require(output_temperature == 0.0_dp, &
      trim(label) // " temperature changed")
    call require(output_sound_speed == 0.0_dp, &
      trim(label) // " sound speed changed")
  end subroutine require_primitive_rejected

  subroutine require_conserved_rejected( &
      active_species, input, temperature_guess, label)
    type(nasa7_species), intent(in) :: active_species(:)
    real(dp), intent(in) :: input(:)
    real(dp), intent(in) :: temperature_guess
    character(len=*), intent(in) :: label

    real(dp), allocatable :: output(:)
    real(dp) :: output_temperature, output_sound_speed
    logical :: local_ok

    allocate(output(reactive_nprim(size(active_species))))
    output = 7.0_dp
    output_temperature = 7.0_dp
    output_sound_speed = 7.0_dp
    call reactive_conserved_to_primitive( &
      active_species, input, temperature_guess, output, output_temperature, &
      output_sound_speed, local_ok)
    call require(.not. local_ok, trim(label) // " accepted")
    call require(all(output == 0.0_dp), trim(label) // " output changed")
    call require(output_temperature == 0.0_dp, &
      trim(label) // " temperature changed")
    call require(output_sound_speed == 0.0_dp, &
      trim(label) // " sound speed changed")
  end subroutine require_conserved_rejected

  subroutine check_2d_timestep_failure_contract( &
      active_species, valid_conserved, valid_temperature)
    type(nasa7_species), intent(in) :: active_species(:)
    real(dp), intent(in) :: valid_conserved(:), valid_temperature

    real(dp), allocatable :: state(:, :, :), temperature_field(:, :)
    real(dp) :: dt, nan_value, positive_inf
    logical :: local_ok
    integer :: i, j

    allocate(state(size(valid_conserved), 2, 2), temperature_field(2, 2))
    do j = 1, 2
      do i = 1, 2
        state(:, i, j) = valid_conserved
        temperature_field(i, j) = valid_temperature
      end do
    end do

    call compute_reactive_cfl_timestep_2d( &
      active_species, state, temperature_field, 2, 2, 0.01_dp, 0.02_dp, &
      0.5_dp, dt, local_ok)
    call require(local_ok .and. ieee_is_finite(dt) .and. dt > 0.0_dp, &
      "valid 2D CFL timestep rejected")

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    positive_inf = ieee_value(0.0_dp, ieee_positive_inf)
    call require_timestep_rejected( &
      active_species, state, temperature_field, nan_value, 0.02_dp, &
      0.5_dp, "2D CFL NaN dx")
    call require_timestep_rejected( &
      active_species, state, temperature_field, 0.01_dp, positive_inf, &
      0.5_dp, "2D CFL infinite dy")
    call require_timestep_rejected( &
      active_species, state, temperature_field, 0.01_dp, 0.02_dp, &
      positive_inf, "2D CFL infinite coefficient")
    call require_timestep_rejected( &
      active_species, state, temperature_field, tiny(1.0_dp), 0.02_dp, &
      0.5_dp, "2D CFL overflowing inverse dx")

  end subroutine check_2d_timestep_failure_contract

  subroutine require_timestep_rejected( &
      active_species, state, temperature_field, dx, dy, cfl, label)
    type(nasa7_species), intent(in) :: active_species(:)
    real(dp), intent(in) :: state(:, :, :), temperature_field(:, :)
    real(dp), intent(in) :: dx, dy, cfl
    character(len=*), intent(in) :: label

    real(dp) :: dt
    logical :: local_ok

    dt = 7.0_dp
    call compute_reactive_cfl_timestep_2d( &
      active_species, state, temperature_field, 2, 2, dx, dy, cfl, dt, &
      local_ok)
    call require(.not. local_ok, trim(label) // " accepted")
    call require(dt == 0.0_dp, trim(label) // " output changed")
  end subroutine require_timestep_rejected

end program test_reactive_general_eos
