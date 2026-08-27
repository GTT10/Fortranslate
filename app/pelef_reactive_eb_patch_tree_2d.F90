program pelef_reactive_eb_patch_tree_2d
  use precision_mod, only: dp
  use constants_mod, only: pelef_version
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use transport_database_mod, only: &
    gas_transport_species, load_h2o2_elementary_transport, &
    load_h2o2_full_transport
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use h2o2_elementary_mechanism_mod, only: &
    load_h2o2_elementary_mechanism
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use simulation_config_reactive_eb_amr_2d_mod, only: &
    reactive_eb_amr_2d_config, read_reactive_eb_amr_2d_configuration
  use amr_eb_patch_tree_reactive_2d_mod, only: &
    reactive_amr_eb_patch_tree_2d, write_reactive_amr_eb_patch_tree_2d_csv
  use reactive_eb_amr_2d_driver_mod, only: &
    simulate_reactive_amr_eb_patch_tree_2d
  implicit none

  type(reactive_eb_amr_2d_config) :: config
  type(reactive_amr_eb_patch_tree_2d) :: solution
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  real(dp), allocatable :: initial_integrals(:), final_integrals(:)
  integer, allocatable :: chemistry_level_advances(:)
  integer, allocatable :: transport_level_advances(:)
  integer, allocatable :: hydro_level_advances(:)
  real(dp) :: base_density, conservation_error, minimum_dt
  real(dp) :: minimum_transport_theta, time
  character(len=1024) :: input_path, message
  integer :: cumulative_tagged_cells, level, regrid_evaluations
  integer :: regrids, steps
  logical :: ok

  if (command_argument_count() /= 1) then
    write(*, '(a)') "Usage: pelef_reactive_eb_patch_tree_2d <input.nml>"
    error stop 2
  end if
  call get_command_argument(1, input_path)
  call read_reactive_eb_amr_2d_configuration( &
    trim(input_path), config, ok, message)
  if (.not. ok) then
    write(*, '(a)') trim(message)
    error stop 2
  end if
  if (config%three_level_enabled .or. config%multipatch_enabled) &
    error stop "Patch-tree application excludes legacy fixed-depth modes"

  select case (trim(config%eb%flow%chemistry_model))
  case ("elementary")
    call load_h2o2_elementary_thermo(species, ok)
    if (.not. ok) error stop "Failed to load elementary thermodynamics"
    call load_h2o2_elementary_mechanism(reactions, ok)
    if (.not. ok) error stop "Failed to load elementary mechanism"
    call load_h2o2_elementary_transport(transport, ok)
    if (.not. ok) error stop "Failed to load elementary transport"
  case ("full_h2o2")
    call load_h2o2_full_thermo(species, ok)
    if (.not. ok) error stop "Failed to load full H2/O2 thermodynamics"
    call load_h2o2_full_mechanism(reactions, ok)
    if (.not. ok) error stop "Failed to load full H2/O2 mechanism"
    call load_h2o2_full_transport(transport, ok)
    if (.not. ok) error stop "Failed to load full H2/O2 transport"
  case default
    error stop "Unknown chemistry model"
  end select

  call simulate_reactive_amr_eb_patch_tree_2d( &
    species, reactions, transport, config, solution, time, steps, regrids, &
    initial_integrals, final_integrals, minimum_dt, base_density, ok, &
    message, minimum_transport_theta, chemistry_level_advances, &
    transport_level_advances, hydro_level_advances, &
    regrid_evaluations, cumulative_tagged_cells)
  if (.not. ok) then
    write(*, '(a,1x,a)') "Patch-tree failure stage:", trim(message)
    error stop "Reactive EB patch-tree 2D simulation failed"
  end if
  call write_reactive_amr_eb_patch_tree_2d_csv( &
    config%eb%flow%output_file, species, solution, time, ok)
  if (.not. ok) error stop "Reactive EB patch-tree composite output failed"

  conservation_error = maxval(abs(final_integrals - initial_integrals) / &
    max(1.0_dp, abs(initial_integrals)))
  write(*, '(a)') "PeleF " // pelef_version // " reactive EB patch-tree 2D"
  write(*, '(a,i0)') "Levels: ", solution%level_count()
  do level = 1, solution%level_count()
    write(*, '(a,i0,a,i0)') &
      "Level ", level - 1, " patches: ", &
      solution%levels(level)%patch_count()
  end do
  write(*, '(a,i0)') "Completed root steps: ", steps
  write(*, '(a,i0)') "Completed regrids: ", regrids
  write(*, '(a,i0)') "Regrid evaluations: ", regrid_evaluations
  write(*, '(a,i0)') "Cumulative tagged cells: ", cumulative_tagged_cells
  write(*, '(a,l2)') "Stopped after checkpoint: ", &
    time < config%eb%flow%final_time
  write(*, '(a,es24.16)') "Final time: ", time
  write(*, '(a,es24.16)') "Minimum accepted root dt: ", minimum_dt
  write(*, '(a,es24.16)') "Base density: ", base_density
  write(*, '(a,es24.16)') "Minimum transport limiter theta: ", &
    minimum_transport_theta
  write(*, '(a,es24.16)') "Maximum composite conservation error: ", &
    conservation_error
  write(*, '(a,*(1x,i0))') "Chemistry level advances:", &
    chemistry_level_advances
  write(*, '(a,*(1x,i0))') "Transport level advances:", &
    transport_level_advances
  write(*, '(a,*(1x,i0))') "Hydro level advances:", hydro_level_advances
  write(*, '(a,1x,a)') "Composite output:", &
    trim(config%eb%flow%output_file)
  if (len_trim(config%checkpoint_file) > 0) &
    write(*, '(a,1x,a)') "Checkpoint:", trim(config%checkpoint_file)

end program pelef_reactive_eb_patch_tree_2d
