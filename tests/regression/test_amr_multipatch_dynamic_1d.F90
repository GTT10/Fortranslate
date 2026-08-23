program test_amr_multipatch_dynamic_1d
  use precision_mod, only: dp
  use state_indices_mod, only: irho, imx, iet
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use simulation_config_reactive_1d_mod, only: reactive_1d_config
  use reactive_1d_mod, only: reactive_nvar
  use amr_multipatch_reactive_1d_mod, only: &
    amr_multipatch_reactive_solution_1d, &
    initialize_tagged_multipatch_reactive_1d, &
    regrid_multipatch_reactive_1d, multipatch_reactive_integrals_1d
  implicit none

  real(dp), parameter :: conservation_tolerance = 8.0e-12_dp
  type(nasa7_species), allocatable :: species(:)
  type(reactive_1d_config) :: config
  type(amr_multipatch_reactive_solution_1d) :: solution
  real(dp), allocatable :: before(:), after(:)
  real(dp) :: error
  logical :: changed, ok

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "dynamic multipatch thermodynamics load")
  call configure_case(config)
  call initialize_tagged_multipatch_reactive_1d( &
    species, config, solution, ok)
  call require(ok .and. solution%patch_count() == 0, &
    "uniform field starts with an empty patch set")
  call require(solution%regrid_evaluations == 1 .and. &
    solution%regrids == 0, "initial empty tag accounting")
  allocate(before(reactive_nvar(size(species))))
  allocate(after(reactive_nvar(size(species))))

  call set_velocity_boxes(solution, 8, 12, 34, 38, 20.0_dp)
  call multipatch_reactive_integrals_1d(solution, before, ok)
  call require(ok, "pre-creation composite integral")
  call regrid_multipatch_reactive_1d( &
    species, config, solution, changed, ok)
  call require(ok .and. changed, "two-patch creation changes the hierarchy")
  call require(solution%patch_count() == 2, &
    "two separated velocity features create two patches")
  call check_integral(before, solution, "patch creation conservation")

  call set_velocity_boxes(solution, 10, 14, 32, 36, 20.0_dp)
  call multipatch_reactive_integrals_1d(solution, before, ok)
  call require(ok, "pre-movement composite integral")
  call regrid_multipatch_reactive_1d( &
    species, config, solution, changed, ok)
  call require(ok .and. changed, "shifted features move the patch set")
  call require(solution%patch_count() == 2, &
    "moving features retain two separated patches")
  call require(solution%overlap_cells_transferred > 0, &
    "patch movement transfers aligned fine overlap")
  call check_integral(before, solution, "patch movement conservation")

  call set_velocity_boxes(solution, 1, 0, 1, 0, 0.0_dp)
  call multipatch_reactive_integrals_1d(solution, before, ok)
  call require(ok, "pre-removal composite integral")
  call regrid_multipatch_reactive_1d( &
    species, config, solution, changed, ok)
  call require(ok .and. changed, "uniform field removes the patch set")
  call require(solution%patch_count() == 0, &
    "empty tags remove every fine patch")
  call check_integral(before, solution, "patch removal conservation")
  call require(solution%regrid_evaluations == 4 .and. &
    solution%regrids == 3, "dynamic patch-set accounting")

  write(*, '(a,1x,es16.8)') &
    "Dynamic multipatch maximum conservation error:", error
  write(*, '(a)') "test_amr_multipatch_dynamic_1d: PASS"

contains

  subroutine configure_case(local_config)
    type(reactive_1d_config), intent(out) :: local_config

    local_config = reactive_1d_config()
    local_config%nx = 48
    local_config%x_lower = 0.0_dp
    local_config%x_upper = 0.012_dp
    local_config%problem = "uniform_reactor"
    local_config%boundary_condition = "periodic"
    local_config%chemistry_enabled = .false.
    local_config%transport_enabled = .false.
    local_config%initial_temperature = 1200.0_dp
    local_config%initial_pressure = 101325.0_dp
    local_config%initial_velocity = 0.0_dp
    local_config%amr_enabled = .true.
    local_config%amr_multipatch_enabled = .true.
    local_config%amr_reconstruction = "plm"
    local_config%amr_refinement_ratio = 2
    local_config%amr_max_levels = 2
    local_config%amr_regrid_interval = 1
    local_config%amr_tag_component = imx
    local_config%amr_buffer_cells = 1
    local_config%amr_minimum_patch_cells = 4
    local_config%amr_maximum_patch_gap_cells = 4
    local_config%amr_relative_gradient_threshold = 0.18_dp
    local_config%amr_absolute_gradient_threshold = 1.0e-6_dp
    local_config%amr_scale_floor = 1.0e-3_dp
  end subroutine configure_case

  subroutine set_velocity_boxes( &
      local_solution, first_lower, first_upper, second_lower, second_upper, &
      amplitude)
    type(amr_multipatch_reactive_solution_1d), intent(inout) :: local_solution
    integer, intent(in) :: first_lower, first_upper
    integer, intent(in) :: second_lower, second_upper
    real(dp), intent(in) :: amplitude

    real(dp) :: x, patch_lower
    integer :: patch, cell, fine_cells, parent_cell

    do cell = 1, local_solution%hierarchy%coarse_cells
      call set_cell_velocity(local_solution%coarse(:, cell), cell, &
        first_lower, first_upper, second_lower, second_upper, amplitude)
    end do
    do patch = 1, local_solution%patch_count()
      patch_lower = local_solution%hierarchy%x_lower + real( &
        local_solution%hierarchy%patches(patch)%fine_coarse_lower - 1, dp) * &
        local_solution%hierarchy%coarse_dx
      fine_cells = &
        local_solution%hierarchy%patches(patch)%fine%cell_count()
      do cell = 1, fine_cells
        x = patch_lower + (real(cell, dp) - 0.5_dp) * &
          local_solution%hierarchy%fine_dx
        parent_cell = min(local_solution%hierarchy%coarse_cells, max(1, &
          int((x - local_solution%hierarchy%x_lower) / &
            local_solution%hierarchy%coarse_dx) + 1))
        call set_cell_velocity(local_solution%patches(patch)%state(:, cell), &
          parent_cell, first_lower, first_upper, second_lower, second_upper, &
          amplitude)
      end do
    end do
  end subroutine set_velocity_boxes

  subroutine set_cell_velocity( &
      state, parent_cell, first_lower, first_upper, second_lower, &
      second_upper, amplitude)
    real(dp), intent(inout) :: state(:)
    integer, intent(in) :: parent_cell, first_lower, first_upper
    integer, intent(in) :: second_lower, second_upper
    real(dp), intent(in) :: amplitude

    real(dp) :: rho, old_momentum, velocity

    rho = state(irho)
    old_momentum = state(imx)
    velocity = 0.0_dp
    if ((parent_cell >= first_lower .and. parent_cell <= first_upper) .or. &
        (parent_cell >= second_lower .and. parent_cell <= second_upper)) &
      velocity = amplitude
    state(imx) = rho * velocity
    state(iet) = state(iet) - 0.5_dp * old_momentum**2 / rho + &
      0.5_dp * state(imx)**2 / rho
  end subroutine set_cell_velocity

  subroutine check_integral(reference, local_solution, label)
    real(dp), intent(in) :: reference(:)
    type(amr_multipatch_reactive_solution_1d), intent(in) :: local_solution
    character(len=*), intent(in) :: label

    call multipatch_reactive_integrals_1d(local_solution, after, ok)
    call require(ok, trim(label) // " evaluation")
    error = maxval(abs(after - reference) / max(1.0_dp, abs(reference)))
    call require(error < conservation_tolerance, label)
  end subroutine check_integral

  subroutine require(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label

    if (.not. condition) then
      write(*, '(a,1x,a)') "FAIL:", trim(label)
      error stop 1
    end if
  end subroutine require

end program test_amr_multipatch_dynamic_1d
