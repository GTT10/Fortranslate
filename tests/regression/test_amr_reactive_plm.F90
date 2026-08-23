program test_amr_reactive_plm
  use precision_mod, only: dp
  use state_indices_mod, only: irho
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: &
    load_h2o2_elementary_mechanism
  use simulation_config_reactive_1d_mod, only: reactive_1d_config
  use reactive_1d_mod, only: initialize_reactive_1d
  use amr_reactive_1d_mod, only: &
    amr_reactive_solution_1d, simulate_amr_reactive_1d
  implicit none

  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(reactive_1d_config) :: config, reference_config
  type(amr_reactive_solution_1d) :: pcm_solution, plm_solution
  real(dp), allocatable :: reference_coarse(:, :), reference_fine(:, :)
  real(dp), allocatable :: reference_temperature(:)
  real(dp) :: pcm_initial(5), pcm_final(5), plm_initial(5), plm_final(5)
  real(dp) :: pcm_error, plm_error, dx
  logical :: ok

  call load_h2o2_elementary_thermo(species, ok)
  if (.not. ok) error stop "Failed to load PLM thermodynamics"
  call load_h2o2_elementary_mechanism(reactions, ok)
  if (.not. ok) error stop "Failed to load PLM mechanism"

  config = reactive_1d_config()
  config%nx = 32
  config%x_lower = 0.0_dp
  config%x_upper = 0.012_dp
  config%final_time = 2.0e-6_dp
  config%cfl = 0.25_dp
  config%maximum_steps = 5000
  config%problem = "reactive_hotspot"
  config%reconstruction = "pcm"
  config%riemann_solver = "rusanov"
  config%limiter = "mc"
  config%boundary_condition = "periodic"
  config%chemistry_enabled = .false.
  config%transport_enabled = .false.
  config%initial_temperature = 1200.0_dp
  config%initial_pressure = 101325.0_dp
  config%initial_velocity = 100.0_dp
  config%hotspot_temperature_rise = 250.0_dp
  config%hotspot_center = 0.006_dp
  config%hotspot_width = 0.0012_dp
  config%amr_enabled = .true.
  config%amr_refinement_ratio = 2
  config%amr_regrid_interval = 1
  config%amr_tag_component = irho
  config%amr_buffer_cells = 2
  config%amr_minimum_patch_cells = 8
  config%amr_relative_gradient_threshold = 0.003_dp
  config%amr_absolute_gradient_threshold = 0.0_dp
  config%amr_scale_floor = 1.0e-12_dp

  config%amr_reconstruction = "pcm"
  call simulate_amr_reactive_1d( &
    species, reactions, config, pcm_solution, pcm_initial, pcm_final, ok)
  if (.not. ok) error stop "AMR PCM comparison run failed"
  config%amr_reconstruction = "plm"
  call simulate_amr_reactive_1d( &
    species, reactions, config, plm_solution, plm_initial, plm_final, ok)
  if (.not. ok) error stop "AMR PLM comparison run failed"
  if (.not. pcm_solution%fine_active() .or. &
      .not. plm_solution%fine_active()) then
    error stop "AMR PLM comparison lost the fine patch"
  end if

  reference_config = config
  reference_config%amr_enabled = .false.
  reference_config%hotspot_center = config%hotspot_center + &
    config%initial_velocity * config%final_time
  reference_config%nx = config%nx
  call initialize_reactive_1d( &
    species, reference_config, reference_coarse, reference_temperature, &
    dx, ok)
  if (.not. ok) error stop "Coarse exact contact initialization failed"
  reference_config%nx = config%nx * config%amr_refinement_ratio
  call initialize_reactive_1d( &
    species, reference_config, reference_fine, reference_temperature, dx, ok)
  if (.not. ok) error stop "Fine exact contact initialization failed"

  call composite_density_error( &
    pcm_solution, reference_coarse, reference_fine, pcm_error)
  call composite_density_error( &
    plm_solution, reference_coarse, reference_fine, plm_error)
  if (plm_error >= 0.85_dp * pcm_error) then
    write(*, '(a,2(1x,es16.8))') "PCM/PLM errors:", pcm_error, plm_error
    error stop "AMR PLM did not improve the advected contact"
  end if
  if (maxval(abs(pcm_final - pcm_initial) / &
      max(1.0_dp, abs(pcm_initial))) > 2.0e-10_dp) then
    error stop "AMR PCM comparison lost conservation"
  end if
  if (maxval(abs(plm_final - plm_initial) / &
      max(1.0_dp, abs(plm_initial))) > 2.0e-10_dp) then
    error stop "AMR PLM comparison lost conservation"
  end if

  write(*, '(a,1x,es16.8)') "AMR PCM density L1:", pcm_error
  write(*, '(a,1x,es16.8)') "AMR PLM density L1:", plm_error
  write(*, '(a)') "test_amr_reactive_plm: PASS"

contains

  subroutine composite_density_error( &
      solution, exact_coarse, exact_fine, error)
    type(amr_reactive_solution_1d), intent(in) :: solution
    real(dp), intent(in) :: exact_coarse(:, 0:), exact_fine(:, 0:)
    real(dp), intent(out) :: error
    real(dp) :: numerator, denominator, weight
    integer :: cell, global_fine, fine_cells

    numerator = 0.0_dp
    denominator = 0.0_dp
    do cell = 1, solution%hierarchy%fine_coarse_lower - 1
      weight = solution%coarse_dx
      numerator = numerator + weight * abs( &
        solution%coarse(irho, cell) - exact_coarse(irho, cell))
      denominator = denominator + weight * abs(exact_coarse(irho, cell))
    end do
    fine_cells = solution%hierarchy%fine%cell_count()
    do cell = 1, fine_cells
      global_fine = solution%hierarchy%fine%lower + cell - 1
      weight = solution%hierarchy%fine_dx
      numerator = numerator + weight * abs( &
        solution%fine(irho, cell) - exact_fine(irho, global_fine))
      denominator = denominator + weight * abs(exact_fine(irho, global_fine))
    end do
    do cell = solution%hierarchy%fine_coarse_upper + 1, config%nx
      weight = solution%coarse_dx
      numerator = numerator + weight * abs( &
        solution%coarse(irho, cell) - exact_coarse(irho, cell))
      denominator = denominator + weight * abs(exact_coarse(irho, cell))
    end do
    error = numerator / denominator
  end subroutine composite_density_error

end program test_amr_reactive_plm
