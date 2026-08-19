program test_reactive_hotspot_reference
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use simulation_config_reactive_1d_mod, only: reactive_1d_config
  use reactive_1d_mod, only: &
    simulate_reactive_1d, reactive_conserved_to_primitive, &
    reactive_nprim, reactive_mass_fraction_component
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  real(dp), allocatable :: reference_state(:, :), reference_temperature(:)
  real(dp) :: plm_error(2), pcm_error(2)
  real(dp) :: dx, time, initial_integrals(5), final_integrals(5)
  logical :: ok
  integer :: steps

  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok) error stop "Failed to load hotspot-reference thermodynamics"
  call load_h2o2_elementary_mechanism(reactions, ok)
  if (.not. ok) error stop "Failed to load hotspot-reference mechanism"

  call run_case(128, "characteristic_plm", reference_state, &
    reference_temperature, dx, time, steps, initial_integrals, &
    final_integrals, ok)
  if (.not. ok) error stop "Reactive hotspot reference run failed"

  call compare_grid(32, "characteristic_plm", reference_state, &
    reference_temperature, plm_error(1))
  call compare_grid(64, "characteristic_plm", reference_state, &
    reference_temperature, plm_error(2))
  call compare_grid(32, "pcm", reference_state, reference_temperature, &
    pcm_error(1))
  call compare_grid(64, "pcm", reference_state, reference_temperature, &
    pcm_error(2))

  write(*, '(a,2(1x,es16.8))') "hotspot PLM normalized L1:", plm_error
  write(*, '(a,2(1x,es16.8))') "hotspot PCM normalized L1:", pcm_error
  if (plm_error(2) >= 0.70_dp * plm_error(1)) then
    error stop "Reactive hotspot PLM did not improve under refinement"
  end if
  if (plm_error(1) >= 0.75_dp * pcm_error(1) .or. &
      plm_error(2) >= 0.75_dp * pcm_error(2)) then
    error stop "Reactive characteristic PLM did not beat PCM"
  end if
  write(*, '(a)') "test_reactive_hotspot_reference: PASS"

contains

  subroutine set_config(nx, reconstruction, config)
    integer, intent(in) :: nx
    character(len=*), intent(in) :: reconstruction
    type(reactive_1d_config), intent(out) :: config

    config = reactive_1d_config()
    config%nx = nx
    config%x_lower = 0.0_dp
    config%x_upper = 0.012_dp
    config%final_time = 8.0e-6_dp
    config%cfl = 0.35_dp
    config%maximum_steps = 50000
    config%problem = "reactive_hotspot"
    config%reconstruction = trim(reconstruction)
    config%limiter = "mc"
    config%boundary_condition = "periodic"
    config%chemistry_enabled = .true.
    config%chemistry_relative_tolerance = 2.0e-7_dp
    config%chemistry_absolute_tolerance = 1.0e-12_dp
    config%initial_temperature = 1200.0_dp
    config%initial_pressure = 101325.0_dp
    config%initial_velocity = 0.0_dp
    config%hotspot_temperature_rise = 250.0_dp
    config%hotspot_center = 0.006_dp
    config%hotspot_width = 0.0012_dp
    config%x_h2 = 0.29570_dp
    config%x_h = 1.0e-5_dp
    config%x_o = 1.0e-5_dp
    config%x_o2 = 0.14784_dp
    config%x_oh = 1.0e-5_dp
    config%x_h2o = 0.0_dp
    config%x_n2 = 0.55643_dp
  end subroutine set_config

  subroutine run_case(nx, reconstruction, state, temperature, dx, time, &
      steps, initial_integrals, final_integrals, run_ok)
    integer, intent(in) :: nx
    character(len=*), intent(in) :: reconstruction
    real(dp), allocatable, intent(out) :: state(:, :), temperature(:)
    real(dp), intent(out) :: dx, time
    integer, intent(out) :: steps
    real(dp), intent(out) :: initial_integrals(5), final_integrals(5)
    logical, intent(out) :: run_ok
    type(reactive_1d_config) :: config

    call set_config(nx, reconstruction, config)
    call simulate_reactive_1d(species, reactions, config, state, temperature, &
      dx, time, steps, initial_integrals, final_integrals, run_ok)
  end subroutine run_case

  subroutine compare_grid(nx, reconstruction, reference, reference_t, error)
    integer, intent(in) :: nx
    character(len=*), intent(in) :: reconstruction
    real(dp), intent(in) :: reference(:, 0:), reference_t(0:)
    real(dp), intent(out) :: error
    real(dp), allocatable :: state(:, :), temperature(:)
    real(dp), allocatable :: q(:), qref(:)
    real(dp) :: dx_local, time_local
    real(dp) :: initial_local(5), final_local(5)
    real(dp) :: coarse_values(5), reference_values(5)
    real(dp) :: scales(5), total_error
    logical :: run_ok, local_ok
    integer :: steps_local, ratio, cell, fine, first, last, k

    call run_case(nx, reconstruction, state, temperature, dx_local, &
      time_local, steps_local, initial_local, final_local, run_ok)
    if (.not. run_ok) error stop "Reactive hotspot comparison run failed"
    ratio = 128 / nx
    if (ratio * nx /= 128) error stop "Invalid hotspot reference ratio"
    allocate(q(reactive_nprim(size(species))))
    allocate(qref(reactive_nprim(size(species))))
    scales = 0.0_dp
    total_error = 0.0_dp
    do cell = 1, nx
      call state_values(state(:, cell), temperature(cell), q, coarse_values, &
        local_ok)
      if (.not. local_ok) error stop "Invalid coarse hotspot state"
      reference_values = 0.0_dp
      first = (cell - 1) * ratio + 1
      last = cell * ratio
      do fine = first, last
        call state_values(reference(:, fine), reference_t(fine), qref, &
          scales, local_ok)
        if (.not. local_ok) error stop "Invalid hotspot reference state"
        reference_values = reference_values + scales
      end do
      reference_values = reference_values / real(ratio, dp)
      scales(1) = max(abs(reference_values(1)), 1.0e-6_dp)
      scales(2) = max(abs(reference_values(2)), 1.0_dp)
      scales(3) = max(abs(reference_values(3)), 1.0e4_dp)
      scales(4) = max(abs(reference_values(4)), 300.0_dp)
      scales(5) = max(abs(reference_values(5)), 1.0e-6_dp)
      do k = 1, 5
        total_error = total_error + &
          abs(coarse_values(k) - reference_values(k)) / scales(k)
      end do
    end do
    error = total_error / real(5 * nx, dp)
    if (maxval(abs(final_local - initial_local) / &
        max(1.0_dp, abs(initial_local))) > 1.0e-11_dp) then
      error stop "Reactive hotspot comparison lost conservation"
    end if
  end subroutine compare_grid

  subroutine state_values(conserved, temperature_guess, primitive, values, ok)
    real(dp), intent(in) :: conserved(:), temperature_guess
    real(dp), intent(out) :: primitive(:), values(5)
    logical, intent(out) :: ok
    real(dp) :: local_temperature, sound_speed

    call reactive_conserved_to_primitive(species, conserved, &
      temperature_guess, primitive, local_temperature, sound_speed, ok)
    if (.not. ok) return
    values(1) = primitive(1)
    values(2) = primitive(2)
    values(3) = primitive(5)
    values(4) = local_temperature
    values(5) = primitive(reactive_mass_fraction_component(5))
  end subroutine state_values

end program test_reactive_hotspot_reference
