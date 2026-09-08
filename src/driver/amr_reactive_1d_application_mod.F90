module amr_reactive_1d_application_mod
  use precision_mod, only: dp
  use constants_mod, only: pelef_version
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use simulation_config_reactive_1d_mod, only: reactive_1d_config
  use amr_reactive_1d_mod, only: &
    amr_reactive_solution_1d, simulate_amr_reactive_1d, &
    write_amr_reactive_1d_csv
  use amr_multilevel_reactive_1d_mod, only: &
    amr_multilevel_reactive_solution_1d, &
    simulate_multilevel_reactive_1d, write_multilevel_reactive_1d_csv
  use amr_multipatch_reactive_1d_mod, only: &
    amr_multipatch_reactive_solution_1d, &
    simulate_multipatch_reactive_1d, write_multipatch_reactive_1d_csv
  implicit none
  private

  public :: run_amr_reactive_1d_application

contains

  subroutine run_amr_reactive_1d_application( &
      input_path, application_label, bundle_sha256, config, species, &
      reactions, transport, base_mole_fractions, chemistry_integrator)
    character(len=*), intent(in) :: input_path, application_label
    character(len=*), intent(in) :: bundle_sha256
    type(reactive_1d_config), intent(in) :: config
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: base_mole_fractions(:)
    character(len=*), intent(in), optional :: chemistry_integrator

    type(amr_reactive_solution_1d) :: solution
    type(amr_multilevel_reactive_solution_1d) :: multilevel_solution
    type(amr_multipatch_reactive_solution_1d) :: multipatch_solution
    real(dp) :: initial_integrals(5), final_integrals(5)
    real(dp) :: conservation_error(5)
    logical :: ok, multilevel_run, multipatch_run
    integer :: patch

    if (.not. config%amr_enabled .or. size(species) < 2 .or. &
        size(reactions) < 1 .or. size(transport) /= size(species) .or. &
        size(base_mole_fractions) /= size(species)) then
      error stop "AMR reactive 1D application data are incompatible"
    end if

    multipatch_run = config%amr_multipatch_enabled
    multilevel_run = .not. multipatch_run .and. &
      (config%amr_max_levels > 2 .or. &
      trim(config%amr_reconstruction) == "ppm" .or. &
      trim(config%amr_reconstruction) == "characteristic_ppm")
    if (multipatch_run) then
      call simulate_multipatch_reactive_1d( &
        species, reactions, config, multipatch_solution, initial_integrals, &
        final_integrals, ok, transport, base_mole_fractions, &
        chemistry_integrator=chemistry_integrator)
      if (.not. ok) error stop "Multipatch AMR reactive simulation failed"
      call write_multipatch_reactive_1d_csv( &
        config%output_file, species, multipatch_solution, ok)
      if (.not. ok) error stop "Multipatch AMR reactive output failed"
    else if (multilevel_run) then
      call simulate_multilevel_reactive_1d( &
        species, reactions, config, multilevel_solution, initial_integrals, &
        final_integrals, ok, transport, base_mole_fractions, &
        chemistry_integrator=chemistry_integrator)
      if (.not. ok) error stop "Multilevel AMR reactive simulation failed"
      call write_multilevel_reactive_1d_csv( &
        config%output_file, species, multilevel_solution, ok)
      if (.not. ok) error stop "Multilevel AMR reactive output failed"
    else
      call simulate_amr_reactive_1d( &
        species, reactions, config, solution, initial_integrals, &
        final_integrals, ok, transport, base_mole_fractions, &
        chemistry_integrator=chemistry_integrator)
      if (.not. ok) error stop "AMR reactive 1D simulation failed"
      call write_amr_reactive_1d_csv( &
        config%output_file, species, solution, ok)
      if (.not. ok) error stop "AMR reactive 1D output failed"
    end if

    conservation_error = abs(final_integrals - initial_integrals) / &
      max(1.0_dp, abs(initial_integrals))
    write(*, '(a)') &
      "PeleF " // pelef_version // " " // trim(application_label)
    if (len_trim(bundle_sha256) > 0) then
      write(*, '(a,1x,a)') "Bundle SHA-256:", trim(bundle_sha256)
      write(*, '(a,i0)') "Species: ", size(species)
      write(*, '(a,i0)') "Reactions: ", size(reactions)
      if (present(chemistry_integrator)) then
        write(*, '(a,1x,a)') "Chemistry integrator:", &
          trim(chemistry_integrator)
      end if
    end if
    write(*, '(a,1x,a)') "Input:", trim(input_path)
    write(*, '(a,1x,a)') "Problem:", trim(config%problem)
    write(*, '(a,i0)') "Coarse cells: ", config%nx
    write(*, '(a,i0)') "Refinement ratio: ", config%amr_refinement_ratio
    write(*, '(a,i0)') "Maximum AMR levels: ", config%amr_max_levels
    write(*, '(a,1x,a)') "AMR reconstruction:", &
      trim(config%amr_reconstruction)
    write(*, '(a,l2)') "AMR hybrid WENO: ", config%amr_hybrid_weno
    if (config%amr_hybrid_weno) then
      write(*, '(a,i0)') "AMR WENO scheme: ", config%amr_weno_scheme
    end if
    write(*, '(a,l2)') "AMR molecular transport: ", &
      config%transport_enabled
    write(*, '(a,l2)') "Dynamic multipatch AMR: ", multipatch_run
    if (multipatch_run) then
      write(*, '(a,i0)') "Active fine patches: ", &
        multipatch_solution%patch_count()
      do patch = 1, multipatch_solution%patch_count()
        write(*, '(a,i0,a,i0,a,i0)') "Patch ", patch, &
          " coarse-cell bounds: ", &
          multipatch_solution%hierarchy%patches(patch)%fine_coarse_lower, &
          ":", &
          multipatch_solution%hierarchy%patches(patch)%fine_coarse_upper
      end do
      write(*, '(a,i0)') "Fine cells: ", &
        multipatch_solution%hierarchy%fine_cell_count()
      write(*, '(a,i0)') "Completed coarse steps: ", &
        multipatch_solution%steps
      write(*, '(a,i0)') "Regrid evaluations: ", &
        multipatch_solution%regrid_evaluations
      write(*, '(a,i0)') "Patch-set changes: ", multipatch_solution%regrids
      write(*, '(a,i0)') "Fine overlap cells transferred: ", &
        multipatch_solution%overlap_cells_transferred
      write(*, '(a,es24.16)') "Final time: ", multipatch_solution%time
    else if (multilevel_run) then
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
      write(*, '(a,i0)') "Regrid evaluations: ", &
        solution%regrid_evaluations
      write(*, '(a,i0)') "Hierarchy changes: ", solution%regrids
      write(*, '(a,es24.16)') "Final time: ", solution%time
    end if
    write(*, '(a,es24.16)') "Maximum conservation error: ", &
      maxval(conservation_error)
    write(*, '(a,1x,a)') "Output:", trim(config%output_file)
  end subroutine run_amr_reactive_1d_application

end module amr_reactive_1d_application_mod
