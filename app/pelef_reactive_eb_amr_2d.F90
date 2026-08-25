program pelef_reactive_eb_amr_2d
  use precision_mod, only: dp
  use constants_mod, only: pelef_version
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use eb_geometry_2d_mod, only: &
    eb_geometry_2d, eb_covered_cell, eb_cut_cell, eb_regular_cell
  use amr_eb_hierarchy_2d_mod, only: amr_eb_patch_2d
  use simulation_config_reactive_eb_2d_mod, only: reactive_eb_2d_config
  use simulation_config_reactive_eb_amr_2d_mod, only: &
    reactive_eb_amr_2d_config, read_reactive_eb_amr_2d_configuration
  use reactive_eb_2d_driver_mod, only: write_reactive_eb_2d_csv
  use reactive_eb_amr_2d_driver_mod, only: simulate_reactive_eb_amr_2d
  implicit none

  type(reactive_eb_amr_2d_config) :: config
  type(reactive_eb_2d_config) :: fine_output_config
  type(eb_geometry_2d) :: coarse_geometry, fine_geometry
  type(amr_eb_patch_2d) :: patch
  type(nasa7_species), allocatable :: species(:)
  real(dp), allocatable :: coarse_state(:, :, :), coarse_temperature(:, :)
  real(dp), allocatable :: fine_state(:, :, :), fine_temperature(:, :)
  real(dp), allocatable :: initial_integrals(:), final_integrals(:)
  real(dp) :: time, minimum_dt, base_density, conservation_error
  character(len=1024) :: input_path, message
  logical :: fine_active, ok
  integer :: regrids, steps

  if (command_argument_count() /= 1) then
    write(*, '(a)') "Usage: pelef_reactive_eb_amr_2d <input.nml>"
    error stop 2
  end if
  call get_command_argument(1, input_path)
  call read_reactive_eb_amr_2d_configuration( &
    trim(input_path), config, ok, message)
  if (.not. ok) then
    write(*, '(a)') trim(message)
    error stop 2
  end if

  select case (trim(config%eb%flow%chemistry_model))
  case ("elementary")
    call load_h2o2_elementary_thermo(species, ok)
  case ("full_h2o2")
    call load_h2o2_full_thermo(species, ok)
  case default
    error stop "Unknown chemistry model"
  end select
  if (.not. ok) error stop "Failed to load thermodynamics"

  call simulate_reactive_eb_amr_2d( &
    species, config, coarse_state, coarse_temperature, coarse_geometry, &
    fine_state, fine_temperature, fine_geometry, patch, fine_active, time, &
    steps, regrids, initial_integrals, final_integrals, minimum_dt, &
    base_density, ok)
  if (.not. ok) error stop "Reactive EB AMR 2D simulation failed"
  call write_reactive_eb_2d_csv( &
    config%eb%flow%output_file, species, config%eb, coarse_state, &
    coarse_temperature, coarse_geometry, time, ok)
  if (.not. ok) error stop "Reactive EB AMR coarse output failed"
  if (fine_active) then
    fine_output_config = config%eb
    fine_output_config%flow%nx = fine_geometry%nx
    fine_output_config%flow%ny = fine_geometry%ny
    fine_output_config%flow%x_lower = fine_geometry%x_lower
    fine_output_config%flow%x_upper = fine_geometry%x_upper
    fine_output_config%flow%y_lower = fine_geometry%y_lower
    fine_output_config%flow%y_upper = fine_geometry%y_upper
    fine_output_config%flow%output_file = trim(config%fine_output_file)
    call write_reactive_eb_2d_csv( &
      config%fine_output_file, species, fine_output_config, fine_state, &
      fine_temperature, fine_geometry, time, ok)
    if (.not. ok) error stop "Reactive EB AMR fine output failed"
  end if

  conservation_error = maxval(abs(final_integrals - initial_integrals) / &
    max(1.0_dp, abs(initial_integrals)))
  write(*, '(a)') "PeleF " // pelef_version // " reactive EB AMR 2D"
  write(*, '(a,i0,a,i0)') &
    "Coarse grid: ", coarse_geometry%nx, " x ", coarse_geometry%ny
  if (fine_active) then
    write(*, '(a,i0,a,i0)') &
      "Fine grid: ", fine_geometry%nx, " x ", fine_geometry%ny
    write(*, '(a,4(i0,1x))') "Coarse patch bounds: ", &
      patch%coarse_i_lower, patch%coarse_i_upper, &
      patch%coarse_j_lower, patch%coarse_j_upper
    write(*, '(a,i0)') "Refinement ratio: ", patch%refinement_ratio
  else
    write(*, '(a)') "Fine grid: inactive"
  end if
  write(*, '(a,i0)') "Coarse regular cells: ", &
    count(coarse_geometry%cell_type == eb_regular_cell)
  write(*, '(a,i0)') "Coarse cut cells: ", &
    count(coarse_geometry%cell_type == eb_cut_cell)
  write(*, '(a,i0)') "Coarse covered cells: ", &
    count(coarse_geometry%cell_type == eb_covered_cell)
  if (fine_active) then
    write(*, '(a,i0)') "Fine regular cells: ", &
      count(fine_geometry%cell_type == eb_regular_cell)
    write(*, '(a,i0)') "Fine cut cells: ", &
      count(fine_geometry%cell_type == eb_cut_cell)
    write(*, '(a,i0)') "Fine covered cells: ", &
      count(fine_geometry%cell_type == eb_covered_cell)
  end if
  write(*, '(a,i0)') "Completed coarse steps: ", steps
  write(*, '(a,i0)') "Completed regrids: ", regrids
  write(*, '(a,es24.16)') "Final time: ", time
  write(*, '(a,es24.16)') "Minimum accepted coarse dt: ", minimum_dt
  write(*, '(a,es24.16)') "Maximum composite conservation error: ", &
    conservation_error
  write(*, '(a,1x,a)') "Coarse output:", &
    trim(config%eb%flow%output_file)
  if (fine_active) then
    write(*, '(a,1x,a)') "Fine output:", trim(config%fine_output_file)
  else
    write(*, '(a)') "Fine output: inactive"
  end if
end program pelef_reactive_eb_amr_2d
