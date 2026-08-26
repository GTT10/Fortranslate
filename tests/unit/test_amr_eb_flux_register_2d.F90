program test_amr_eb_flux_register_2d
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_primitive_to_conserved
  use eb_geometry_2d_mod, only: &
    eb_geometry_2d, eb_covered_cell, eb_cut_cell, build_eb_geometry_2d
  use amr_eb_hierarchy_2d_mod, only: &
    amr_eb_patch_2d, build_amr_eb_patch_2d, composite_eb_integral_2d
  use amr_eb_flux_register_2d_mod, only: &
    amr_eb_flux_register_2d, initialize_amr_eb_flux_register_2d, &
    accumulate_coarse_eb_fluxes_2d, accumulate_fine_eb_fluxes_2d, &
    reflux_eb_state_patch_2d, &
    reflux_reactive_eb_state_patch_support_2d, &
    reflux_reactive_eb_state_patch_2d
  implicit none

  integer, parameter :: coarse_nx = 8, coarse_ny = 8
  integer, parameter :: coarse_i_lower = 2, coarse_i_upper = 6
  integer, parameter :: coarse_j_lower = 2, coarse_j_upper = 6
  integer, parameter :: ratio = 2
  integer, parameter :: fine_nx = &
    (coarse_i_upper - coarse_i_lower + 1) * ratio
  integer, parameter :: fine_ny = &
    (coarse_j_upper - coarse_j_lower + 1) * ratio
  type(eb_geometry_2d) :: coarse_geometry, fine_geometry
  type(amr_eb_patch_2d) :: patch
  type(amr_eb_flux_register_2d) :: register, support_register
  type(nasa7_species), allocatable :: species(:)
  real(dp) :: coarse_level_set(0:coarse_nx, 0:coarse_ny)
  real(dp) :: fine_level_set(0:fine_nx, 0:fine_ny)
  real(dp), allocatable :: coarse_x_flux(:, :, :), coarse_y_flux(:, :, :)
  real(dp), allocatable :: fine_x_flux(:, :, :), fine_y_flux(:, :, :)
  real(dp), allocatable :: coarse_state(:, :, :), fine_state(:, :, :)
  real(dp), allocatable :: refluxed_coarse(:, :, :), refluxed_fine(:, :, :)
  real(dp), allocatable :: raw_correction(:, :, :), composite_integral(:)
  real(dp), allocatable :: reactive_coarse_x(:, :, :), reactive_coarse_y(:, :, :)
  real(dp), allocatable :: reactive_fine_x(:, :, :), reactive_fine_y(:, :, :)
  real(dp), allocatable :: reactive_coarse(:, :, :), reactive_fine(:, :, :)
  real(dp), allocatable :: reactive_refluxed_coarse(:, :, :)
  real(dp), allocatable :: reactive_refluxed_fine(:, :, :)
  real(dp), allocatable :: support_refluxed_coarse(:, :, :)
  real(dp), allocatable :: support_refluxed_fine(:, :, :)
  real(dp), allocatable :: coarse_temperature(:, :), fine_temperature(:, :)
  real(dp), allocatable :: refluxed_coarse_temperature(:, :)
  real(dp), allocatable :: refluxed_fine_temperature(:, :)
  real(dp), allocatable :: support_coarse_temperature(:, :)
  real(dp), allocatable :: support_fine_temperature(:, :)
  real(dp), allocatable :: primitive(:), state_cell(:), mass_fractions(:)
  real(dp) :: mole_fractions(7), x, y, fine_x_lower, fine_x_upper
  real(dp) :: fine_y_lower, fine_y_upper, dt, raw_integral
  real(dp) :: temperature_cell, sound_speed, state_scale
  logical :: ok
  integer :: i, j, nvar

  do j = 0, coarse_ny
    y = real(j, dp) / real(coarse_ny, dp)
    do i = 0, coarse_nx
      x = real(i, dp) / real(coarse_nx, dp)
      coarse_level_set(i, j) = x + y - 0.78_dp
    end do
  end do
  call build_eb_geometry_2d( &
    coarse_level_set, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    coarse_geometry, ok)
  call require(ok, "coarse EB geometry")

  fine_x_lower = real(coarse_i_lower - 1, dp) / real(coarse_nx, dp)
  fine_x_upper = real(coarse_i_upper, dp) / real(coarse_nx, dp)
  fine_y_lower = real(coarse_j_lower - 1, dp) / real(coarse_ny, dp)
  fine_y_upper = real(coarse_j_upper, dp) / real(coarse_ny, dp)
  do j = 0, fine_ny
    y = fine_y_lower + real(j, dp) * &
      (fine_y_upper - fine_y_lower) / real(fine_ny, dp)
    do i = 0, fine_nx
      x = fine_x_lower + real(i, dp) * &
        (fine_x_upper - fine_x_lower) / real(fine_nx, dp)
      fine_level_set(i, j) = x + y - 0.78_dp
    end do
  end do
  call build_eb_geometry_2d( &
    fine_level_set, fine_x_lower, fine_x_upper, &
    fine_y_lower, fine_y_upper, fine_geometry, ok)
  call require(ok, "fine EB geometry")
  call build_amr_eb_patch_2d( &
    coarse_geometry, fine_geometry, coarse_i_lower, coarse_i_upper, &
    coarse_j_lower, coarse_j_upper, ratio, patch, ok)
  call require(ok, "EB AMR patch")
  call require(coarse_geometry%cell_type(1, 6) == eb_cut_cell, &
    "qualified coarse/fine cut cell")

  allocate(coarse_x_flux(2, 0:coarse_nx, coarse_ny))
  allocate(coarse_y_flux(2, coarse_nx, 0:coarse_ny))
  allocate(fine_x_flux(2, 0:fine_nx, fine_ny))
  allocate(fine_y_flux(2, fine_nx, 0:fine_ny))
  call initialize_amr_eb_flux_register_2d( &
    coarse_geometry, fine_geometry, patch, 2, register, ok)
  call require(ok .and. register%is_valid( &
    coarse_geometry, fine_geometry, patch), "flux-register initialization")
  call require( &
    register%correction_i_lower == coarse_i_lower - 1 .and. &
    register%correction_i_upper == coarse_i_upper + 1 .and. &
    register%correction_j_lower == coarse_j_lower - 1 .and. &
    register%correction_j_upper == coarse_j_upper + 1 .and. &
    lbound(register%correction, 2) == coarse_i_lower - 1 .and. &
    ubound(register%correction, 2) == coarse_i_upper + 1 .and. &
    lbound(register%correction, 3) == coarse_j_lower - 1 .and. &
    ubound(register%correction, 3) == coarse_j_upper + 1 .and. &
    size(register%correction, 2) * size(register%correction, 3) < &
      coarse_nx * coarse_ny, &
    "compact flux-register correction support")

  coarse_x_flux(1, :, :) = 1.2_dp
  coarse_x_flux(2, :, :) = -0.7_dp
  coarse_y_flux(1, :, :) = 0.4_dp
  coarse_y_flux(2, :, :) = 0.9_dp
  fine_x_flux(1, :, :) = 1.2_dp
  fine_x_flux(2, :, :) = -0.7_dp
  fine_y_flux(1, :, :) = 0.4_dp
  fine_y_flux(2, :, :) = 0.9_dp
  dt = 0.04_dp
  call accumulate_coarse_eb_fluxes_2d( &
    register, coarse_geometry, fine_geometry, patch, &
    coarse_x_flux, coarse_y_flux, dt, ok)
  call require(ok, "coarse EB flux accumulation")
  call accumulate_fine_eb_fluxes_2d( &
    register, coarse_geometry, fine_geometry, patch, &
    fine_x_flux, fine_y_flux, 0.5_dp * dt, ok)
  call require(ok, "first fine EB flux accumulation")
  call accumulate_fine_eb_fluxes_2d( &
    register, coarse_geometry, fine_geometry, patch, &
    fine_x_flux, fine_y_flux, 0.5_dp * dt, ok)
  call require(ok, "second fine EB flux accumulation")
  call require(maxval(abs(register%correction)) <= 3.0e-14_dp, &
    "aperture-matched subcycled flux cancellation")

  call register%reset()
  coarse_x_flux = 0.0_dp
  coarse_y_flux = 0.0_dp
  fine_x_flux = 0.0_dp
  fine_y_flux = 0.0_dp
  fine_x_flux(1, 0, 9:10) = 1.0_dp
  call accumulate_fine_eb_fluxes_2d( &
    register, coarse_geometry, fine_geometry, patch, &
    fine_x_flux, fine_y_flux, 0.02_dp, ok)
  call require(ok .and. abs(register%correction(1, 1, 6)) > 0.0_dp, &
    "cut-interface flux mismatch accumulation")
  allocate(raw_correction, source=register%correction)
  raw_integral = sum(coarse_geometry%volume_fraction( &
    register%correction_i_lower:register%correction_i_upper, &
    register%correction_j_lower:register%correction_j_upper) * &
    raw_correction(1, :, :)) * coarse_geometry%dx * coarse_geometry%dy

  allocate(coarse_state(2, coarse_nx, coarse_ny))
  allocate(fine_state(2, fine_nx, fine_ny))
  allocate(refluxed_coarse(2, coarse_nx, coarse_ny))
  allocate(refluxed_fine(2, fine_nx, fine_ny))
  allocate(composite_integral(2))
  coarse_state = 0.0_dp
  fine_state = 0.0_dp
  call reflux_eb_state_patch_2d( &
    coarse_state, coarse_geometry, fine_state, fine_geometry, patch, &
    register, refluxed_coarse, refluxed_fine, ok)
  call require(ok, "EB cut-cell re-reflux")
  call composite_eb_integral_2d( &
    refluxed_coarse, coarse_geometry, refluxed_fine, fine_geometry, &
    patch, composite_integral, ok)
  call require(ok, "reflux composite integral")
  call assert_close(composite_integral(1), raw_integral, 8.0e-13_dp, &
    "re-reflux preserves raw extensive correction")
  call assert_close(refluxed_coarse(1, 1, 6), &
    coarse_geometry%volume_fraction(1, 6) * raw_correction(1, 1, 6), &
    4.0e-14_dp, "cut-cell self correction is volume scaled")
  call require(maxval(abs(refluxed_fine)) > 0.0_dp, &
    "fine-covered re-reflux recipient")
  call require(maxval(abs(register%correction)) == 0.0_dp, &
    "successful re-reflux resets register")

  call initialize_amr_eb_flux_register_2d( &
    coarse_geometry, fine_geometry, patch, 2, register, ok)
  call require(ok, "transaction register initialization")
  fine_x_flux = 0.0_dp
  fine_x_flux(1, 0, 9:10) = 1.0_dp
  fine_x_flux(2, 0, 9:10) = -2.0_dp
  call accumulate_fine_eb_fluxes_2d( &
    register, coarse_geometry, fine_geometry, patch, &
    fine_x_flux, fine_y_flux, 0.02_dp, ok)
  call require(ok, "transaction baseline accumulation")
  deallocate(raw_correction)
  allocate(raw_correction, source=register%correction)
  fine_x_flux(1, 0, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
  call accumulate_fine_eb_fluxes_2d( &
    register, coarse_geometry, fine_geometry, patch, &
    fine_x_flux, fine_y_flux, 0.02_dp, ok)
  call require(.not. ok .and. &
    maxval(abs(register%correction - raw_correction)) == 0.0_dp, &
    "nonfinite accumulation transaction")

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "thermodynamic database load")
  nvar = reactive_nvar(size(species))
  allocate(primitive(reactive_nprim(size(species))))
  allocate(mass_fractions(size(species)), state_cell(nvar))
  mole_fractions = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
    1.0e-5_dp, 0.0_dp, 0.55643_dp]
  call mass_fractions_from_mole_fractions( &
    species, mole_fractions, mass_fractions, ok)
  call require(ok, "composition conversion")
  primitive(1:5) = [0.31_dp, 4.0_dp, -2.0_dp, 0.0_dp, 135000.0_dp]
  do i = 1, size(species)
    primitive(reactive_mass_fraction_component(i)) = mass_fractions(i)
  end do
  call reactive_primitive_to_conserved( &
    species, primitive, state_cell, temperature_cell, sound_speed, ok)
  call require(ok, "reference reactive state")
  allocate(reactive_coarse_x(nvar, 0:coarse_nx, coarse_ny))
  allocate(reactive_coarse_y(nvar, coarse_nx, 0:coarse_ny))
  allocate(reactive_fine_x(nvar, 0:fine_nx, fine_ny))
  allocate(reactive_fine_y(nvar, fine_nx, 0:fine_ny))
  allocate(reactive_coarse(nvar, coarse_nx, coarse_ny))
  allocate(reactive_fine(nvar, fine_nx, fine_ny))
  allocate(reactive_refluxed_coarse(nvar, coarse_nx, coarse_ny))
  allocate(reactive_refluxed_fine(nvar, fine_nx, fine_ny))
  allocate(coarse_temperature(coarse_nx, coarse_ny))
  allocate(fine_temperature(fine_nx, fine_ny))
  allocate(refluxed_coarse_temperature(coarse_nx, coarse_ny))
  allocate(refluxed_fine_temperature(fine_nx, fine_ny))
  do j = 1, coarse_ny
    do i = 1, coarse_nx
      reactive_coarse(:, i, j) = state_cell
      coarse_temperature(i, j) = temperature_cell
    end do
  end do
  do j = 1, fine_ny
    do i = 1, fine_nx
      reactive_fine(:, i, j) = state_cell
      fine_temperature(i, j) = temperature_cell
    end do
  end do
  reactive_coarse_x = 0.0_dp
  reactive_coarse_y = 0.0_dp
  reactive_fine_x = 0.0_dp
  reactive_fine_y = 0.0_dp
  reactive_fine_x(:, 0, 9:10) = &
    spread(1.0e-3_dp * state_cell, 2, 2)
  call initialize_amr_eb_flux_register_2d( &
    coarse_geometry, fine_geometry, patch, nvar, register, ok)
  call require(ok, "reactive register initialization")
  call accumulate_fine_eb_fluxes_2d( &
    register, coarse_geometry, fine_geometry, patch, &
    reactive_fine_x, reactive_fine_y, 1.0e-3_dp, ok)
  call require(ok, "reactive fine-flux accumulation")
  support_register = register
  allocate(support_refluxed_coarse, mold=reactive_coarse)
  allocate(support_refluxed_fine, mold=reactive_fine)
  allocate(support_coarse_temperature, mold=coarse_temperature)
  allocate(support_fine_temperature, mold=fine_temperature)
  call reflux_reactive_eb_state_patch_support_2d( &
    species, 1, 1, reactive_coarse, coarse_temperature, coarse_geometry, &
    reactive_fine, fine_temperature, fine_geometry, patch, &
    support_register, support_refluxed_coarse, support_coarse_temperature, &
    support_refluxed_fine, support_fine_temperature, ok)
  call require(ok .and. maxval(abs(support_register%correction)) == 0.0_dp, &
    "reactive support re-reflux transaction")
  call reflux_reactive_eb_state_patch_2d( &
    species, reactive_coarse, coarse_temperature, coarse_geometry, &
    reactive_fine, fine_temperature, fine_geometry, patch, register, &
    reactive_refluxed_coarse, refluxed_coarse_temperature, &
    reactive_refluxed_fine, refluxed_fine_temperature, ok)
  call require(ok, "reactive EB re-reflux transaction")
  call require( &
    maxval(abs(support_refluxed_coarse - &
      reactive_refluxed_coarse)) == 0.0_dp .and. &
    maxval(abs(support_refluxed_fine - &
      reactive_refluxed_fine)) == 0.0_dp .and. &
    maxval(abs(support_coarse_temperature - &
      refluxed_coarse_temperature)) == 0.0_dp .and. &
    maxval(abs(support_fine_temperature - &
      refluxed_fine_temperature)) == 0.0_dp, &
    "reactive support/full re-reflux parity")
  state_scale = max(1.0_dp, maxval(abs(state_cell)))
  call require(maxval(abs(reactive_refluxed_coarse - reactive_coarse)) > &
    1.0e-15_dp * state_scale, "reactive re-reflux changes active state")
  call require(maxval(abs(refluxed_coarse_temperature - coarse_temperature)) <= &
    2.0e-8_dp .and. &
    maxval(abs(refluxed_fine_temperature - fine_temperature)) <= 2.0e-8_dp, &
    "scaled-state temperature recovery")
  call assert_covered_unchanged( &
    reactive_coarse, reactive_refluxed_coarse, coarse_geometry, &
    "coarse covered state preservation")
  call assert_covered_unchanged( &
    reactive_fine, reactive_refluxed_fine, fine_geometry, &
    "fine covered state preservation")
  call require(maxval(abs(register%correction)) == 0.0_dp, &
    "reactive success resets register")

  reactive_fine_x = 0.0_dp
  reactive_fine_x(:, 0, 9:10) = spread(1.0e6_dp * state_cell, 2, 2)
  call initialize_amr_eb_flux_register_2d( &
    coarse_geometry, fine_geometry, patch, nvar, register, ok)
  call require(ok, "reactive failure register initialization")
  call accumulate_fine_eb_fluxes_2d( &
    register, coarse_geometry, fine_geometry, patch, &
    reactive_fine_x, reactive_fine_y, 1.0_dp, ok)
  call require(ok .and. maxval(abs(register%correction)) > 0.0_dp, &
    "nonphysical register accumulation")
  call reflux_reactive_eb_state_patch_2d( &
    species, reactive_coarse, coarse_temperature, coarse_geometry, &
    reactive_fine, fine_temperature, fine_geometry, patch, register, &
    reactive_refluxed_coarse, refluxed_coarse_temperature, &
    reactive_refluxed_fine, refluxed_fine_temperature, ok)
  call require(.not. ok .and. &
    maxval(abs(reactive_refluxed_coarse - reactive_coarse)) == 0.0_dp .and. &
    maxval(abs(reactive_refluxed_fine - reactive_fine)) == 0.0_dp .and. &
    maxval(abs(refluxed_coarse_temperature - coarse_temperature)) == 0.0_dp .and. &
    maxval(abs(refluxed_fine_temperature - fine_temperature)) == 0.0_dp .and. &
    maxval(abs(register%correction)) > 0.0_dp, &
    "nonphysical reactive re-reflux rollback")

  write(*, '(a)') "test_amr_eb_flux_register_2d: PASS"

contains

  subroutine assert_covered_unchanged( &
      original, candidate, geometry, message)
    real(dp), intent(in) :: original(:, :, :), candidate(:, :, :)
    type(eb_geometry_2d), intent(in) :: geometry
    character(len=*), intent(in) :: message
    integer :: cell_i, cell_j

    do cell_j = 1, geometry%ny
      do cell_i = 1, geometry%nx
        if (geometry%cell_type(cell_i, cell_j) /= eb_covered_cell) cycle
        call require(maxval(abs(candidate(:, cell_i, cell_j) - &
          original(:, cell_i, cell_j))) == 0.0_dp, message)
      end do
    end do
  end subroutine assert_covered_unchanged

  subroutine assert_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: message

    call require(abs(actual - expected) <= &
      tolerance * max(1.0_dp, abs(expected)), message)
  end subroutine assert_close

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_amr_eb_flux_register_2d
