module reactive_eb_amr_2d_application_mod
  use precision_mod, only: dp
  use constants_mod, only: pelef_version
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use eb_geometry_2d_mod, only: &
    eb_geometry_2d, eb_covered_cell, eb_cut_cell, eb_regular_cell
  use amr_eb_hierarchy_2d_mod, only: amr_eb_patch_2d
  use amr_eb_regrid_2d_mod, only: reactive_eb_patch_set_2d
  use simulation_config_reactive_eb_2d_mod, only: reactive_eb_2d_config
  use simulation_config_reactive_eb_amr_2d_mod, only: &
    reactive_eb_amr_2d_config
  use reactive_eb_2d_driver_mod, only: write_reactive_eb_2d_csv
  use reactive_eb_amr_2d_driver_mod, only: &
    simulate_reactive_eb_amr_2d, simulate_reactive_eb_amr_patch_set_2d, &
    simulate_three_level_reactive_eb_amr_2d
  implicit none
  private

  public :: run_reactive_eb_amr_2d_application

contains

  subroutine run_reactive_eb_amr_2d_application( &
      application_label, bundle_sha256, config, species, reactions, &
      transport, base_mole_fractions, chemistry_integrator)
    character(len=*), intent(in) :: application_label, bundle_sha256
    type(reactive_eb_amr_2d_config), intent(in) :: config
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in), optional :: base_mole_fractions(:)
    character(len=*), intent(in), optional :: chemistry_integrator

    type(reactive_eb_2d_config) :: fine_output_config
    type(eb_geometry_2d) :: coarse_geometry, fine_geometry
    type(eb_geometry_2d) :: level_two_geometry
    type(amr_eb_patch_2d) :: patch, level_two_patch
    type(reactive_eb_patch_set_2d) :: patch_set
    real(dp), allocatable :: coarse_state(:, :, :), coarse_temperature(:, :)
    real(dp), allocatable :: fine_state(:, :, :), fine_temperature(:, :)
    real(dp), allocatable :: level_two_state(:, :, :)
    real(dp), allocatable :: level_two_temperature(:, :)
    real(dp), allocatable :: initial_integrals(:), final_integrals(:)
    real(dp) :: time, minimum_dt, base_density, conservation_error
    real(dp) :: minimum_transport_theta
    character(len=1024) :: message, patch_output_file
    logical :: fine_active, ok, selected_context
    integer :: child, regrids, steps

    if (size(species) < 2 .or. size(reactions) < 1 .or. &
        size(transport) /= size(species)) then
      error stop "Reactive EB AMR 2D application data have incompatible dimensions"
    end if
    if (present(base_mole_fractions)) then
      if (size(base_mole_fractions) /= size(species)) then
        error stop "Reactive EB AMR 2D composition has incompatible dimensions"
      end if
    end if
    selected_context = len_trim(bundle_sha256) > 0
    if (selected_context .neqv. &
        (present(base_mole_fractions) .and. &
         present(chemistry_integrator))) then
      error stop "Reactive EB AMR 2D selected context is incomplete"
    end if
    if ((present(base_mole_fractions) .neqv. &
         present(chemistry_integrator)) .or. &
        (.not. selected_context .and. &
         (present(base_mole_fractions) .or. &
          present(chemistry_integrator)))) then
      error stop "Reactive EB AMR 2D selected context is inconsistent"
    end if
    message = ""

    if (config%three_level_enabled) then
      if (selected_context) then
        call simulate_three_level_reactive_eb_amr_2d( &
          species, reactions, config, coarse_state, coarse_temperature, &
          coarse_geometry, fine_state, fine_temperature, fine_geometry, &
          patch, level_two_state, level_two_temperature, level_two_geometry, &
          level_two_patch, time, steps, regrids, initial_integrals, &
          final_integrals, minimum_dt, base_density, ok, message, transport, &
          minimum_transport_theta, base_mole_fractions, &
          chemistry_integrator, bundle_sha256)
      else
        call simulate_three_level_reactive_eb_amr_2d( &
          species, reactions, config, coarse_state, coarse_temperature, &
          coarse_geometry, fine_state, fine_temperature, fine_geometry, &
          patch, level_two_state, level_two_temperature, level_two_geometry, &
          level_two_patch, time, steps, regrids, initial_integrals, &
          final_integrals, minimum_dt, base_density, ok, message, transport, &
          minimum_transport_theta)
      end if
      fine_active = .true.
    else if (config%multipatch_enabled) then
      if (selected_context) then
        call simulate_reactive_eb_amr_patch_set_2d( &
          species, reactions, config, coarse_state, coarse_temperature, &
          coarse_geometry, patch_set, time, steps, regrids, &
          initial_integrals, final_integrals, minimum_dt, base_density, ok, &
          message, transport, minimum_transport_theta, base_mole_fractions, &
          chemistry_integrator, bundle_sha256)
      else
        call simulate_reactive_eb_amr_patch_set_2d( &
          species, reactions, config, coarse_state, coarse_temperature, &
          coarse_geometry, patch_set, time, steps, regrids, &
          initial_integrals, final_integrals, minimum_dt, base_density, ok, &
          message, transport, minimum_transport_theta)
      end if
      fine_active = .false.
    else
      if (selected_context) then
        call simulate_reactive_eb_amr_2d( &
          species, reactions, config, coarse_state, coarse_temperature, &
          coarse_geometry, fine_state, fine_temperature, fine_geometry, &
          patch, fine_active, time, steps, regrids, initial_integrals, &
          final_integrals, minimum_dt, base_density, ok, transport, &
          minimum_transport_theta, base_mole_fractions, &
          chemistry_integrator, bundle_sha256, message)
      else
        call simulate_reactive_eb_amr_2d( &
          species, reactions, config, coarse_state, coarse_temperature, &
          coarse_geometry, fine_state, fine_temperature, fine_geometry, &
          patch, fine_active, time, steps, regrids, initial_integrals, &
          final_integrals, minimum_dt, base_density, ok, transport, &
          minimum_transport_theta, failure_context=message)
      end if
    end if
    if (.not. ok) then
      if (config%three_level_enabled) &
        write(*, '(a,1x,a)') "Three-level failure stage:", trim(message)
      if (config%multipatch_enabled) &
        write(*, '(a,1x,a)') "Multipatch failure stage:", trim(message)
      if (.not. config%three_level_enabled .and. &
          .not. config%multipatch_enabled .and. len_trim(message) > 0) &
        write(*, '(a)') trim(message)
      error stop "Reactive EB AMR 2D simulation failed"
    end if

    call write_reactive_eb_2d_csv( &
      config%eb%flow%output_file, species, config%eb, coarse_state, &
      coarse_temperature, coarse_geometry, time, ok)
    if (.not. ok) error stop "Reactive EB AMR coarse output failed"
    if (config%three_level_enabled) then
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
      if (.not. ok) error stop "Reactive EB AMR middle output failed"
    else if (config%multipatch_enabled) then
      do child = 1, patch_set%patch_count()
        call make_patch_output_path( &
          config%fine_output_file, child, patch_output_file)
        fine_output_config = config%eb
        fine_output_config%flow%nx = patch_set%children(child)%geometry%nx
        fine_output_config%flow%ny = patch_set%children(child)%geometry%ny
        fine_output_config%flow%x_lower = &
          patch_set%children(child)%geometry%x_lower
        fine_output_config%flow%x_upper = &
          patch_set%children(child)%geometry%x_upper
        fine_output_config%flow%y_lower = &
          patch_set%children(child)%geometry%y_lower
        fine_output_config%flow%y_upper = &
          patch_set%children(child)%geometry%y_upper
        fine_output_config%flow%output_file = trim(patch_output_file)
        call write_reactive_eb_2d_csv( &
          patch_output_file, species, fine_output_config, &
          patch_set%children(child)%state, &
          patch_set%children(child)%temperature, &
          patch_set%children(child)%geometry, time, ok)
        if (.not. ok) error stop "Reactive EB AMR patch output failed"
      end do
    else if (fine_active) then
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
    if (config%three_level_enabled) then
      fine_output_config = config%eb
      fine_output_config%flow%nx = level_two_geometry%nx
      fine_output_config%flow%ny = level_two_geometry%ny
      fine_output_config%flow%x_lower = level_two_geometry%x_lower
      fine_output_config%flow%x_upper = level_two_geometry%x_upper
      fine_output_config%flow%y_lower = level_two_geometry%y_lower
      fine_output_config%flow%y_upper = level_two_geometry%y_upper
      fine_output_config%flow%output_file = trim(config%level_two_output_file)
      call write_reactive_eb_2d_csv( &
        config%level_two_output_file, species, fine_output_config, &
        level_two_state, level_two_temperature, level_two_geometry, time, ok)
      if (.not. ok) error stop "Reactive EB AMR level-two output failed"
    end if

    conservation_error = maxval(abs(final_integrals - initial_integrals) / &
      max(1.0_dp, abs(initial_integrals)))
    write(*, '(a)') &
      "PeleF " // pelef_version // " " // trim(application_label)
    if (len_trim(bundle_sha256) > 0) then
      write(*, '(a,1x,a)') "Bundle SHA-256:", trim(bundle_sha256)
      if (present(chemistry_integrator)) then
        write(*, '(a,1x,a)') "Chemistry integrator:", &
          trim(chemistry_integrator)
      end if
    end if
    write(*, '(a,i0,a,i0)') &
      "Coarse grid: ", coarse_geometry%nx, " x ", coarse_geometry%ny
    if (config%three_level_enabled) then
      write(*, '(a,i0,a,i0)') &
        "Middle grid: ", fine_geometry%nx, " x ", fine_geometry%ny
      write(*, '(a,4(i0,1x))') "Root patch bounds: ", &
        patch%coarse_i_lower, patch%coarse_i_upper, &
        patch%coarse_j_lower, patch%coarse_j_upper
      write(*, '(a,i0,a,i0)') &
        "Finest grid: ", level_two_geometry%nx, " x ", level_two_geometry%ny
      write(*, '(a,4(i0,1x))') "Middle patch bounds: ", &
        level_two_patch%coarse_i_lower, level_two_patch%coarse_i_upper, &
        level_two_patch%coarse_j_lower, level_two_patch%coarse_j_upper
      write(*, '(a,i0)') "Refinement ratio: ", config%refinement_ratio
    else if (config%multipatch_enabled) then
      write(*, '(a,i0)') "Fine patches: ", patch_set%patch_count()
      do child = 1, patch_set%patch_count()
        write(*, '(a,i0,a,i0,a,i0)') "Fine patch ", child, ": ", &
          patch_set%children(child)%geometry%nx, " x ", &
          patch_set%children(child)%geometry%ny
        write(*, '(a,i0,a,4(i0,1x))') "Fine patch ", child, &
          " coarse bounds: ", &
          patch_set%children(child)%patch%coarse_i_lower, &
          patch_set%children(child)%patch%coarse_i_upper, &
          patch_set%children(child)%patch%coarse_j_lower, &
          patch_set%children(child)%patch%coarse_j_upper
      end do
      write(*, '(a,i0)') "Refinement ratio: ", config%refinement_ratio
    else if (fine_active) then
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
    if (config%multipatch_enabled) then
      do child = 1, patch_set%patch_count()
        write(*, '(a,i0,a,i0)') "Fine patch ", child, " regular cells: ", &
          count(patch_set%children(child)%geometry%cell_type == eb_regular_cell)
        write(*, '(a,i0,a,i0)') "Fine patch ", child, " cut cells: ", &
          count(patch_set%children(child)%geometry%cell_type == eb_cut_cell)
        write(*, '(a,i0,a,i0)') "Fine patch ", child, " covered cells: ", &
          count(patch_set%children(child)%geometry%cell_type == eb_covered_cell)
      end do
    else if (fine_active) then
      write(*, '(a,i0)') "Fine regular cells: ", &
        count(fine_geometry%cell_type == eb_regular_cell)
      write(*, '(a,i0)') "Fine cut cells: ", &
        count(fine_geometry%cell_type == eb_cut_cell)
      write(*, '(a,i0)') "Fine covered cells: ", &
        count(fine_geometry%cell_type == eb_covered_cell)
    end if
    if (config%three_level_enabled) then
      write(*, '(a,i0)') "Finest regular cells: ", &
        count(level_two_geometry%cell_type == eb_regular_cell)
      write(*, '(a,i0)') "Finest cut cells: ", &
        count(level_two_geometry%cell_type == eb_cut_cell)
      write(*, '(a,i0)') "Finest covered cells: ", &
        count(level_two_geometry%cell_type == eb_covered_cell)
    end if
    write(*, '(a,i0)') "Completed coarse steps: ", steps
    write(*, '(a,i0)') "Completed regrids: ", regrids
    write(*, '(a,l2)') "Stopped after checkpoint: ", &
      time < config%eb%flow%final_time
    write(*, '(a,l2)') "Chemistry: ", config%eb%flow%chemistry_enabled
    write(*, '(a,1x,a)') "Chemistry model:", &
      trim(config%eb%flow%chemistry_model)
    write(*, '(a,l2)') "Molecular transport: ", &
      config%eb%flow%transport_enabled
    if (config%eb%flow%transport_enabled) then
      write(*, '(a,l2)') "Viscosity: ", config%eb%flow%viscosity_enabled
      write(*, '(a,l2)') "Thermal conduction: ", &
        config%eb%flow%thermal_conduction_enabled
      write(*, '(a,l2)') "Species diffusion: ", &
        config%eb%flow%species_diffusion_enabled
      write(*, '(a,l2)') "Barodiffusion: ", &
        config%eb%flow%barodiffusion_enabled
      write(*, '(a,es24.16)') "Minimum transport limiter theta: ", &
        minimum_transport_theta
    end if
    write(*, '(a,es24.16)') "Final time: ", time
    write(*, '(a,es24.16)') "Minimum accepted coarse dt: ", minimum_dt
    write(*, '(a,es24.16)') "Maximum composite conservation error: ", &
      conservation_error
    write(*, '(a,1x,a)') "Coarse output:", &
      trim(config%eb%flow%output_file)
    if (config%multipatch_enabled) then
      do child = 1, patch_set%patch_count()
        call make_patch_output_path( &
          config%fine_output_file, child, patch_output_file)
        write(*, '(a,i0,a,1x,a)') &
          "Fine patch ", child, " output:", trim(patch_output_file)
      end do
      if (patch_set%patch_count() == 0) &
        write(*, '(a)') "Fine output: inactive"
    else if (fine_active) then
      write(*, '(a,1x,a)') "Fine output:", trim(config%fine_output_file)
    else
      write(*, '(a)') "Fine output: inactive"
    end if
    if (config%three_level_enabled) &
      write(*, '(a,1x,a)') "Finest output:", &
        trim(config%level_two_output_file)
    if (len_trim(config%checkpoint_file) > 0) &
      write(*, '(a,1x,a)') "Checkpoint:", trim(config%checkpoint_file)
    if (len_trim(config%restart_file) > 0) &
      write(*, '(a,1x,a)') "Restart source:", trim(config%restart_file)

  contains

    subroutine make_patch_output_path(base_path, patch_index, output_path)
      character(len=*), intent(in) :: base_path
      integer, intent(in) :: patch_index
      character(len=*), intent(out) :: output_path

      character(len=16) :: index_text
      integer :: dot

      output_path = ""
      write(index_text, '(i4.4)') patch_index
      dot = scan(trim(base_path), ".", back=.true.)
      if (dot > 1) then
        output_path = trim(base_path(:dot - 1)) // "_patch" // &
          trim(index_text) // trim(base_path(dot:))
      else
        output_path = trim(base_path) // "_patch" // trim(index_text) // ".csv"
      end if
    end subroutine make_patch_output_path

  end subroutine run_reactive_eb_amr_2d_application

end module reactive_eb_amr_2d_application_mod
