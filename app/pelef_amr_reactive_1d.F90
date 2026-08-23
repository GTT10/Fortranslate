program pelef_amr_reactive_1d
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
  use simulation_config_reactive_1d_mod, only: &
    reactive_1d_config, read_reactive_1d_configuration
  use amr_reactive_1d_mod, only: &
    amr_reactive_solution_1d, simulate_amr_reactive_1d, &
    write_amr_reactive_1d_csv
  use amr_multilevel_reactive_1d_mod, only: &
    amr_multilevel_reactive_solution_1d, &
    simulate_multilevel_reactive_1d, write_multilevel_reactive_1d_csv
  implicit none

  type(reactive_1d_config) :: config
  type(amr_reactive_solution_1d) :: solution
  type(amr_multilevel_reactive_solution_1d) :: multilevel_solution
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  real(dp) :: initial_integrals(5), final_integrals(5)
  real(dp) :: conservation_error(5)
  character(len=1024) :: input_path, message
  logical :: ok, multilevel_run

  if (command_argument_count() /= 1) then
    write(*, '(a)') "Usage: pelef_amr_reactive_1d <input.nml>"
    error stop 2
  end if
  call get_command_argument(1, input_path)
  call read_reactive_1d_configuration(trim(input_path), config, ok, message)
  if (.not. ok) then
    write(*, '(a)') trim(message)
    error stop 2
  end if
  if (.not. config%amr_enabled) &
    error stop "AMR reactive application requires amr_enabled"

  select case (trim(config%chemistry_model))
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

  multilevel_run = config%amr_max_levels > 2
  if (multilevel_run) then
    call simulate_multilevel_reactive_1d( &
      species, reactions, config, multilevel_solution, initial_integrals, &
      final_integrals, ok, transport)
    if (.not. ok) error stop "Multilevel AMR reactive simulation failed"
    call write_multilevel_reactive_1d_csv( &
      config%output_file, species, multilevel_solution, ok)
    if (.not. ok) error stop "Multilevel AMR reactive output failed"
  else
    call simulate_amr_reactive_1d( &
      species, reactions, config, solution, initial_integrals, &
      final_integrals, ok, transport)
    if (.not. ok) error stop "AMR reactive 1D simulation failed"
    call write_amr_reactive_1d_csv( &
      config%output_file, species, solution, ok)
    if (.not. ok) error stop "AMR reactive 1D output failed"
  end if

  conservation_error = abs(final_integrals - initial_integrals) / &
    max(1.0_dp, abs(initial_integrals))
  write(*, '(a)') "PeleF " // pelef_version // " AMR reactive 1D"
  write(*, '(a,1x,a)') "Problem:", trim(config%problem)
  write(*, '(a,i0)') "Coarse cells: ", config%nx
  write(*, '(a,i0)') "Refinement ratio: ", config%amr_refinement_ratio
  write(*, '(a,i0)') "Maximum AMR levels: ", config%amr_max_levels
  write(*, '(a,1x,a)') "AMR reconstruction:", &
    trim(config%amr_reconstruction)
  write(*, '(a,l2)') "AMR molecular transport: ", &
    config%transport_enabled
  if (multilevel_run) then
    write(*, '(a,i0)') "Active AMR levels: ", &
      multilevel_solution%level_count()
    write(*, '(a,i0)') "Completed coarse steps: ", &
      multilevel_solution%steps
    write(*, '(a,i0)') "Regrid evaluations: ", &
      multilevel_solution%regrid_evaluations
    write(*, '(a,i0)') "Hierarchy changes: ", multilevel_solution%regrids
    write(*, '(a,i0)') "Fine overlap cells transferred: ", &
      multilevel_solution%overlap_cells_transferred
    write(*, '(a,es24.16)') "Final time: ", multilevel_solution%time
  else
    write(*, '(a,l2)') "Fine level active: ", solution%fine_active()
    if (solution%fine_active()) then
      write(*, '(a,i0,a,i0)') "Fine coarse-cell bounds: ", &
        solution%hierarchy%fine_coarse_lower, ":", &
        solution%hierarchy%fine_coarse_upper
      write(*, '(a,i0)') "Fine cells: ", &
        solution%hierarchy%fine%cell_count()
    end if
    write(*, '(a,i0)') "Completed coarse steps: ", solution%steps
    write(*, '(a,i0)') "Regrid evaluations: ", solution%regrid_evaluations
    write(*, '(a,i0)') "Hierarchy changes: ", solution%regrids
    write(*, '(a,es24.16)') "Final time: ", solution%time
  end if
  write(*, '(a,es24.16)') "Maximum conservation error: ", &
    maxval(conservation_error)
  write(*, '(a,1x,a)') "Output:", trim(config%output_file)
end program pelef_amr_reactive_1d
