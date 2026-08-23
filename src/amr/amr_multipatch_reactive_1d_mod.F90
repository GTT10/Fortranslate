module amr_multipatch_reactive_1d_mod
  use precision_mod, only: dp
  use state_indices_mod, only: irho, imx, imy, imz, iet
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use transport_database_mod, only: gas_transport_species
  use simulation_config_reactive_1d_mod, only: reactive_1d_config
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_cfl_timestep, &
    reactive_transport_timestep, &
    initialize_reactive_1d, advance_reactive_chemistry
  use amr_hierarchy_1d_mod, only: &
    amr_level_field_1d, amr_flux_register_1d, &
    accumulate_coarse_flux_1d, accumulate_fine_flux_1d
  use amr_multipatch_1d_mod, only: &
    amr_patch_set_1d, initialize_patch_set_1d, prolong_patch_set_1d, &
    average_down_patch_set_1d, &
    initialize_patch_flux_registers_1d, synchronize_patch_set_1d, &
    composite_integral_patch_set_1d
  use amr_reactive_1d_mod, only: &
    amr_ppm_ghost_width, advance_amr_level_1d, advance_transport_level_1d, &
    recover_level_temperatures_1d, fill_physical_ghosts_1d, &
    fill_fine_ghosts_1d, fill_fine_wide_ghosts_1d, write_amr_cell
  use amr_regrid_1d_mod, only: &
    amr_tagging_criteria_1d, amr_regrid_plan_collection_1d, &
    tag_gradient_1d, build_regrid_plan_collection_1d, &
    regrid_patch_set_state_1d
  implicit none
  private

  type, public :: amr_multipatch_reactive_patch_1d
    real(dp), allocatable :: state(:, :)
    real(dp), allocatable :: temperature(:)
    real(dp), allocatable :: left_ghost_state(:, :)
    real(dp), allocatable :: right_ghost_state(:, :)
    real(dp), allocatable :: left_ghost_temperature(:)
    real(dp), allocatable :: right_ghost_temperature(:)
  end type amr_multipatch_reactive_patch_1d

  type, public :: amr_multipatch_reactive_solution_1d
    type(amr_patch_set_1d) :: hierarchy
    real(dp), allocatable :: coarse(:, :)
    real(dp), allocatable :: coarse_temperature(:)
    type(amr_multipatch_reactive_patch_1d), allocatable :: patches(:)
    real(dp) :: time = 0.0_dp
    integer :: steps = 0
    integer :: regrid_evaluations = 0
    integer :: regrids = 0
    integer :: overlap_cells_transferred = 0
  contains
    procedure :: patch_count => multipatch_reactive_patch_count
    procedure :: is_valid => multipatch_reactive_is_valid
  end type amr_multipatch_reactive_solution_1d

  public :: initialize_multipatch_reactive_1d
  public :: initialize_tagged_multipatch_reactive_1d
  public :: multipatch_reactive_timestep_1d
  public :: advance_multipatch_reactive_1d
  public :: advance_multipatch_reactive_hydro_1d
  public :: regrid_multipatch_reactive_1d
  public :: simulate_multipatch_reactive_1d
  public :: multipatch_reactive_integrals_1d
  public :: write_multipatch_reactive_1d_csv

contains

  pure integer function multipatch_reactive_patch_count(self) result(count)
    class(amr_multipatch_reactive_solution_1d), intent(in) :: self

    count = 0
    if (allocated(self%patches)) count = size(self%patches)
  end function multipatch_reactive_patch_count

  pure logical function multipatch_reactive_is_valid(self) result(valid)
    class(amr_multipatch_reactive_solution_1d), intent(in) :: self

    integer :: patch, fine_cells, nvar

    valid = self%hierarchy%is_valid() .and. &
      allocated(self%coarse) .and. allocated(self%coarse_temperature) .and. &
      allocated(self%patches) .and. self%time >= 0.0_dp .and. &
      self%steps >= 0 .and. self%regrid_evaluations >= 0 .and. &
      self%regrids >= 0 .and. self%overlap_cells_transferred >= 0
    if (.not. valid) return
    nvar = size(self%coarse, 1)
    valid = nvar >= 1 .and. lbound(self%coarse, 2) == 0 .and. &
      ubound(self%coarse, 2) == self%hierarchy%coarse_cells + 1 .and. &
      lbound(self%coarse_temperature, 1) == 0 .and. &
      ubound(self%coarse_temperature, 1) == &
        self%hierarchy%coarse_cells + 1 .and. &
      size(self%patches) == self%hierarchy%patch_count()
    if (.not. valid) return
    do patch = 1, size(self%patches)
      fine_cells = self%hierarchy%patches(patch)%fine%cell_count()
      valid = allocated(self%patches(patch)%state) .and. &
        allocated(self%patches(patch)%temperature) .and. &
        allocated(self%patches(patch)%left_ghost_state) .and. &
        allocated(self%patches(patch)%right_ghost_state) .and. &
        allocated(self%patches(patch)%left_ghost_temperature) .and. &
        allocated(self%patches(patch)%right_ghost_temperature)
      if (.not. valid) return
      valid = size(self%patches(patch)%state, 1) == nvar .and. &
        lbound(self%patches(patch)%state, 2) == 0 .and. &
        ubound(self%patches(patch)%state, 2) == fine_cells + 1 .and. &
        lbound(self%patches(patch)%temperature, 1) == 0 .and. &
        ubound(self%patches(patch)%temperature, 1) == fine_cells + 1 .and. &
        size(self%patches(patch)%left_ghost_state, 1) == nvar .and. &
        size(self%patches(patch)%right_ghost_state, 1) == nvar .and. &
        size(self%patches(patch)%left_ghost_state, 2) == &
          amr_ppm_ghost_width .and. &
        size(self%patches(patch)%right_ghost_state, 2) == &
          amr_ppm_ghost_width .and. &
        size(self%patches(patch)%left_ghost_temperature) == &
          amr_ppm_ghost_width .and. &
        size(self%patches(patch)%right_ghost_temperature) == &
          amr_ppm_ghost_width
      if (.not. valid) return
    end do
  end function multipatch_reactive_is_valid

  subroutine initialize_multipatch_reactive_1d( &
      species, config, patch_parent_lower, patch_parent_upper, solution, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    integer, intent(in) :: patch_parent_lower(:), patch_parent_upper(:)
    type(amr_multipatch_reactive_solution_1d), intent(out) :: solution
    logical, intent(out) :: ok

    type(amr_level_field_1d), allocatable :: fields(:)
    real(dp) :: coarse_dx
    logical :: local_ok
    integer :: patch, fine_cells, nvar

    ok = .false.
    if (size(species) < 1 .or. .not. config%amr_enabled) return
    call initialize_reactive_1d( &
      species, config, solution%coarse, solution%coarse_temperature, &
      coarse_dx, local_ok)
    if (.not. local_ok) return
    call initialize_patch_set_1d( &
      config%nx, patch_parent_lower, patch_parent_upper, &
      config%amr_refinement_ratio, config%x_lower, config%x_upper, &
      solution%hierarchy, local_ok)
    if (.not. local_ok) return
    call prolong_patch_set_1d( &
      solution%coarse(:, 1:config%nx), solution%hierarchy, fields, local_ok)
    if (.not. local_ok) return

    nvar = reactive_nvar(size(species))
    allocate(solution%patches(solution%hierarchy%patch_count()))
    do patch = 1, size(solution%patches)
      fine_cells = solution%hierarchy%patches(patch)%fine%cell_count()
      allocate(solution%patches(patch)%state(nvar, 0:fine_cells + 1))
      allocate(solution%patches(patch)%temperature(0:fine_cells + 1))
      allocate(solution%patches(patch)%left_ghost_state( &
        nvar, amr_ppm_ghost_width))
      allocate(solution%patches(patch)%right_ghost_state( &
        nvar, amr_ppm_ghost_width))
      allocate(solution%patches(patch)%left_ghost_temperature( &
        amr_ppm_ghost_width))
      allocate(solution%patches(patch)%right_ghost_temperature( &
        amr_ppm_ghost_width))
      solution%patches(patch)%state = 0.0_dp
      solution%patches(patch)%temperature = 0.0_dp
      solution%patches(patch)%left_ghost_state = 0.0_dp
      solution%patches(patch)%right_ghost_state = 0.0_dp
      solution%patches(patch)%left_ghost_temperature = 0.0_dp
      solution%patches(patch)%right_ghost_temperature = 0.0_dp
      solution%patches(patch)%state(:, 1:fine_cells) = &
        fields(patch)%values
      call recover_level_temperatures_1d( &
        species, solution%patches(patch)%state, &
        solution%patches(patch)%temperature, fine_cells, local_ok)
      if (.not. local_ok) return
    end do
    solution%time = 0.0_dp
    solution%steps = 0
    solution%regrid_evaluations = 0
    solution%regrids = 0
    solution%overlap_cells_transferred = 0
    call refresh_multipatch_ghosts(species, config, solution, local_ok)
    if (.not. local_ok) return
    ok = solution%is_valid() .and. coarse_dx > 0.0_dp
  end subroutine initialize_multipatch_reactive_1d

  subroutine initialize_tagged_multipatch_reactive_1d( &
      species, config, solution, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    type(amr_multipatch_reactive_solution_1d), intent(out) :: solution
    logical, intent(out) :: ok

    type(amr_regrid_plan_collection_1d) :: collection
    real(dp), allocatable :: root_state(:, :), root_temperature(:)
    integer, allocatable :: patch_lower(:), patch_upper(:)
    real(dp) :: root_dx
    logical :: local_ok
    integer :: patch

    ok = .false.
    if (.not. valid_tagged_multipatch_configuration(config, &
        reactive_nvar(size(species)))) return
    call initialize_reactive_1d( &
      species, config, root_state, root_temperature, root_dx, local_ok)
    if (.not. local_ok .or. root_dx <= 0.0_dp) return
    call build_tagged_patch_collection( &
      config, root_state(:, 1:config%nx), collection, local_ok)
    if (.not. local_ok) return
    allocate(patch_lower(collection%patch_count()))
    allocate(patch_upper(collection%patch_count()))
    do patch = 1, collection%patch_count()
      patch_lower(patch) = collection%plans(patch)%patch_lower
      patch_upper(patch) = collection%plans(patch)%patch_upper
    end do
    call initialize_multipatch_reactive_1d( &
      species, config, patch_lower, patch_upper, solution, local_ok)
    if (.not. local_ok) return
    solution%regrid_evaluations = 1
    solution%regrids = merge(1, 0, solution%patch_count() > 0)
    ok = solution%is_valid()
  end subroutine initialize_tagged_multipatch_reactive_1d

  subroutine multipatch_reactive_timestep_1d( &
      species, config, solution, dt, ok, transport)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    type(amr_multipatch_reactive_solution_1d), intent(in) :: solution
    real(dp), intent(out) :: dt
    logical, intent(out) :: ok
    type(gas_transport_species), intent(in), optional :: transport(:)

    real(dp) :: local_dt, transport_dt, maximum_diffusivity
    logical :: local_ok
    integer :: patch, fine_cells

    dt = 0.0_dp
    ok = solution%is_valid()
    if (.not. ok) return
    if (config%transport_enabled .and. .not. present(transport)) then
      ok = .false.
      return
    end if
    call reactive_cfl_timestep( &
      species, solution%coarse, solution%coarse_temperature, &
      solution%hierarchy%coarse_cells, solution%hierarchy%coarse_dx, &
      config%cfl, dt, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    if (config%transport_enabled) then
      call reactive_transport_timestep( &
        species, transport, solution%coarse, &
        solution%coarse_temperature, solution%hierarchy%coarse_cells, &
        solution%hierarchy%coarse_dx, config%transport_cfl, &
        config%viscosity_enabled, config%thermal_conduction_enabled, &
        config%species_diffusion_enabled, transport_dt, &
        maximum_diffusivity, local_ok)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
      dt = min(dt, transport_dt)
    end if
    do patch = 1, solution%patch_count()
      fine_cells = solution%hierarchy%patches(patch)%fine%cell_count()
      call reactive_cfl_timestep( &
        species, solution%patches(patch)%state, &
        solution%patches(patch)%temperature, fine_cells, &
        solution%hierarchy%fine_dx, config%cfl, local_dt, local_ok)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
      dt = min(dt, real(solution%hierarchy%refinement_ratio, dp) * local_dt)
      if (config%transport_enabled) then
        call reactive_transport_timestep( &
          species, transport, solution%patches(patch)%state, &
          solution%patches(patch)%temperature, fine_cells, &
          solution%hierarchy%fine_dx, config%transport_cfl, &
          config%viscosity_enabled, config%thermal_conduction_enabled, &
          config%species_diffusion_enabled, transport_dt, &
          maximum_diffusivity, local_ok)
        if (.not. local_ok) then
          ok = .false.
          return
        end if
        dt = min(dt, &
          real(solution%hierarchy%refinement_ratio**2, dp) * transport_dt)
      end if
    end do
    ok = dt > 0.0_dp
  end subroutine multipatch_reactive_timestep_1d

  subroutine advance_multipatch_reactive_1d( &
      species, reactions, config, dt, solution, ok, transport)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(reactive_1d_config), intent(in) :: config
    real(dp), intent(in) :: dt
    type(amr_multipatch_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok
    type(gas_transport_species), intent(in), optional :: transport(:)

    type(amr_multipatch_reactive_solution_1d) :: backup
    logical :: local_ok

    ok = .false.
    if (dt <= 0.0_dp .or. .not. solution%is_valid()) return
    if (config%transport_enabled .and. .not. present(transport)) return
    backup = solution

    if (config%chemistry_enabled) then
      call advance_multipatch_chemistry( &
        species, reactions, config, 0.5_dp * dt, solution, local_ok)
      if (.not. local_ok) then
        solution = backup
        return
      end if
    end if
    if (config%transport_enabled) then
      call advance_multipatch_transport( &
        species, transport, config, 0.5_dp * dt, solution, local_ok)
      if (.not. local_ok) then
        solution = backup
        return
      end if
    end if
    call advance_multipatch_hydro_impl( &
      species, config, dt, solution, local_ok)
    if (.not. local_ok) then
      solution = backup
      return
    end if
    if (config%transport_enabled) then
      call advance_multipatch_transport( &
        species, transport, config, 0.5_dp * dt, solution, local_ok)
      if (.not. local_ok) then
        solution = backup
        return
      end if
    end if
    if (config%chemistry_enabled) then
      call advance_multipatch_chemistry( &
        species, reactions, config, 0.5_dp * dt, solution, local_ok)
      if (.not. local_ok) then
        solution = backup
        return
      end if
    end if
    call refresh_multipatch_ghosts(species, config, solution, local_ok)
    if (.not. local_ok) then
      solution = backup
      return
    end if
    ok = solution%is_valid()
    if (.not. ok) solution = backup
  end subroutine advance_multipatch_reactive_1d

  subroutine advance_multipatch_reactive_hydro_1d( &
      species, config, dt, solution, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    real(dp), intent(in) :: dt
    type(amr_multipatch_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok

    type(amr_multipatch_reactive_solution_1d) :: backup

    backup = solution
    call advance_multipatch_hydro_impl(species, config, dt, solution, ok)
    if (.not. ok) solution = backup
  end subroutine advance_multipatch_reactive_hydro_1d

  subroutine advance_multipatch_hydro_impl( &
      species, config, dt, solution, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    real(dp), intent(in) :: dt
    type(amr_multipatch_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok

    type(amr_flux_register_1d), allocatable :: registers(:)
    type(amr_level_field_1d), allocatable :: fields(:)
    real(dp), allocatable :: coarse_start(:, :), coarse_end(:, :)
    real(dp), allocatable :: coarse_flux(:, :), fine_flux(:, :)
    real(dp) :: child_dt, alpha
    logical :: local_ok
    integer :: nx, nvar, ratio, substep, patch, fine_cells

    ok = solution%is_valid() .and. dt > 0.0_dp
    if (.not. ok) return
    nx = solution%hierarchy%coarse_cells
    nvar = size(solution%coarse, 1)
    ratio = solution%hierarchy%refinement_ratio
    child_dt = dt / real(ratio, dp)
    allocate(coarse_start(nvar, 0:nx + 1))
    allocate(coarse_end(nvar, 0:nx + 1))
    allocate(coarse_flux(nvar, 0:nx))
    coarse_start = solution%coarse
    call advance_amr_level_1d( &
      species, solution%coarse, solution%coarse_temperature, nx, &
      solution%hierarchy%coarse_dx, dt, config%amr_reconstruction, &
      config%limiter, config%riemann_solver, .true., &
      config%boundary_condition, coarse_flux, local_ok, &
      ppm_contact_steepening=config%ppm_contact_steepening, &
      ppm_shock_flattening=config%ppm_shock_flattening, &
      amr_hybrid_weno=config%amr_hybrid_weno, &
      amr_weno_scheme=config%amr_weno_scheme)
    if (.not. local_ok) return
    coarse_end = solution%coarse

    call initialize_patch_flux_registers_1d( &
      solution%hierarchy, nvar, registers, local_ok)
    if (.not. local_ok) return
    do patch = 1, solution%patch_count()
      call accumulate_coarse_flux_1d( &
        registers(patch), &
        coarse_flux(:, &
          solution%hierarchy%patches(patch)%fine_coarse_lower - 1), &
        coarse_flux(:, &
          solution%hierarchy%patches(patch)%fine_coarse_upper), &
        dt, local_ok)
      if (.not. local_ok) return
    end do

    do substep = 1, ratio
      if (trim(config%amr_reconstruction) == "pcm") then
        alpha = real(substep - 1, dp) / real(ratio, dp)
      else
        alpha = (real(substep, dp) - 0.5_dp) / real(ratio, dp)
      end if
      do patch = 1, solution%patch_count()
        fine_cells = solution%hierarchy%patches(patch)%fine%cell_count()
        call fill_fine_ghosts_1d( &
          species, solution%hierarchy%patches(patch), &
          coarse_start, coarse_end, alpha, &
          solution%patches(patch)%state, &
          solution%patches(patch)%temperature, local_ok, &
          config%boundary_condition)
        if (.not. local_ok) return
        if (uses_ppm_reconstruction(config)) then
          call fill_fine_wide_ghosts_1d( &
            species, solution%hierarchy%patches(patch), &
            coarse_start, coarse_end, alpha, &
            solution%patches(patch)%left_ghost_state, &
            solution%patches(patch)%right_ghost_state, &
            solution%patches(patch)%left_ghost_temperature, &
            solution%patches(patch)%right_ghost_temperature, local_ok, &
            solution%patches(patch)%state, &
            solution%patches(patch)%temperature, &
            config%boundary_condition)
          if (.not. local_ok) return
        end if
        allocate(fine_flux(nvar, 0:fine_cells))
        if (uses_ppm_reconstruction(config)) then
          call advance_amr_level_1d( &
            species, solution%patches(patch)%state, &
            solution%patches(patch)%temperature, fine_cells, &
            solution%hierarchy%fine_dx, child_dt, &
            config%amr_reconstruction, config%limiter, &
            config%riemann_solver, .false., "coarse_fine", &
            fine_flux, local_ok, &
            solution%patches(patch)%left_ghost_state, &
            solution%patches(patch)%right_ghost_state, &
            solution%patches(patch)%left_ghost_temperature, &
            solution%patches(patch)%right_ghost_temperature, &
            config%ppm_contact_steepening, &
            config%ppm_shock_flattening, config%amr_hybrid_weno, &
            config%amr_weno_scheme)
        else
          call advance_amr_level_1d( &
            species, solution%patches(patch)%state, &
            solution%patches(patch)%temperature, fine_cells, &
            solution%hierarchy%fine_dx, child_dt, &
            config%amr_reconstruction, config%limiter, &
            config%riemann_solver, .false., "coarse_fine", &
            fine_flux, local_ok)
        end if
        if (.not. local_ok) return
        call accumulate_fine_flux_1d( &
          registers(patch), fine_flux(:, 0), &
          fine_flux(:, fine_cells), child_dt, local_ok)
        if (.not. local_ok) return
        deallocate(fine_flux)
      end do
    end do

    call extract_patch_fields(solution, fields, local_ok)
    if (.not. local_ok) return
    call synchronize_patch_set_1d( &
      solution%coarse(:, 1:nx), fields, solution%hierarchy, &
      registers, local_ok)
    if (.not. local_ok) return
    call recover_level_temperatures_1d( &
      species, solution%coarse, solution%coarse_temperature, nx, local_ok)
    if (.not. local_ok) return
    call refresh_multipatch_ghosts(species, config, solution, local_ok)
    if (.not. local_ok) return
    solution%time = solution%time + dt
    solution%steps = solution%steps + 1
    ok = solution%is_valid()
  end subroutine advance_multipatch_hydro_impl

  subroutine advance_multipatch_transport( &
      species, transport, config, interval, solution, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(reactive_1d_config), intent(in) :: config
    real(dp), intent(in) :: interval
    type(amr_multipatch_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok

    type(amr_flux_register_1d), allocatable :: registers(:)
    type(amr_level_field_1d), allocatable :: fields(:)
    real(dp), allocatable :: coarse_start(:, :), coarse_end(:, :)
    real(dp), allocatable :: coarse_flux(:, :), fine_flux(:, :)
    real(dp) :: child_interval, alpha, boundary_distance
    logical :: local_ok
    integer :: nx, nvar, ratio, subcycles, substep, patch, fine_cells

    ok = solution%is_valid() .and. interval > 0.0_dp
    if (.not. ok) return
    nx = solution%hierarchy%coarse_cells
    nvar = size(solution%coarse, 1)
    ratio = solution%hierarchy%refinement_ratio
    subcycles = ratio * ratio
    child_interval = interval / real(subcycles, dp)
    boundary_distance = 0.5_dp * &
      (solution%hierarchy%coarse_dx + solution%hierarchy%fine_dx)
    allocate(coarse_start(nvar, 0:nx + 1))
    allocate(coarse_end(nvar, 0:nx + 1))
    allocate(coarse_flux(nvar, 0:nx))
    coarse_start = solution%coarse
    call advance_transport_level_1d( &
      species, transport, solution%coarse, solution%coarse_temperature, &
      nx, solution%hierarchy%coarse_dx, interval, &
      solution%hierarchy%coarse_dx, config, .true., &
      config%boundary_condition, coarse_flux, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    coarse_end = solution%coarse

    call initialize_patch_flux_registers_1d( &
      solution%hierarchy, nvar, registers, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    do patch = 1, solution%patch_count()
      call accumulate_coarse_flux_1d( &
        registers(patch), &
        coarse_flux(:, &
          solution%hierarchy%patches(patch)%fine_coarse_lower - 1), &
        coarse_flux(:, &
          solution%hierarchy%patches(patch)%fine_coarse_upper), &
        interval, local_ok)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
    end do

    do substep = 1, subcycles
      alpha = (real(substep, dp) - 0.5_dp) / real(subcycles, dp)
      do patch = 1, solution%patch_count()
        fine_cells = solution%hierarchy%patches(patch)%fine%cell_count()
        call fill_fine_ghosts_1d( &
          species, solution%hierarchy%patches(patch), &
          coarse_start, coarse_end, alpha, &
          solution%patches(patch)%state, &
          solution%patches(patch)%temperature, local_ok, &
          config%boundary_condition)
        if (.not. local_ok) then
          ok = .false.
          return
        end if
        allocate(fine_flux(nvar, 0:fine_cells))
        call advance_transport_level_1d( &
          species, transport, solution%patches(patch)%state, &
          solution%patches(patch)%temperature, fine_cells, &
          solution%hierarchy%fine_dx, child_interval, boundary_distance, &
          config, .false., "coarse_fine", fine_flux, local_ok)
        if (.not. local_ok) then
          ok = .false.
          return
        end if
        call accumulate_fine_flux_1d( &
          registers(patch), fine_flux(:, 0), &
          fine_flux(:, fine_cells), child_interval, local_ok)
        if (.not. local_ok) then
          ok = .false.
          return
        end if
        deallocate(fine_flux)
      end do
    end do

    call extract_patch_fields(solution, fields, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    call synchronize_patch_set_1d( &
      solution%coarse(:, 1:nx), fields, solution%hierarchy, &
      registers, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    call recover_level_temperatures_1d( &
      species, solution%coarse, solution%coarse_temperature, nx, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    call refresh_multipatch_ghosts(species, config, solution, local_ok)
    ok = local_ok
  end subroutine advance_multipatch_transport

  subroutine advance_multipatch_chemistry( &
      species, reactions, config, interval, solution, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(reactive_1d_config), intent(in) :: config
    real(dp), intent(in) :: interval
    type(amr_multipatch_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok

    type(amr_level_field_1d), allocatable :: fields(:)
    logical :: local_ok
    integer :: nx, patch, fine_cells

    ok = solution%is_valid() .and. interval >= 0.0_dp
    if (.not. ok) return
    nx = solution%hierarchy%coarse_cells
    call advance_reactive_chemistry( &
      species, reactions, solution%coarse, solution%coarse_temperature, &
      nx, interval, config%chemistry_relative_tolerance, &
      config%chemistry_absolute_tolerance, config%boundary_condition, &
      local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    do patch = 1, solution%patch_count()
      fine_cells = solution%hierarchy%patches(patch)%fine%cell_count()
      call advance_reactive_chemistry( &
        species, reactions, solution%patches(patch)%state, &
        solution%patches(patch)%temperature, fine_cells, interval, &
        config%chemistry_relative_tolerance, &
        config%chemistry_absolute_tolerance, "outflow", local_ok)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
    end do
    call extract_patch_fields(solution, fields, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    call average_down_patch_set_1d( &
      fields, solution%hierarchy, solution%coarse(:, 1:nx), local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    call recover_level_temperatures_1d( &
      species, solution%coarse, solution%coarse_temperature, nx, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    call refresh_multipatch_ghosts(species, config, solution, local_ok)
    ok = local_ok
  end subroutine advance_multipatch_chemistry

  subroutine regrid_multipatch_reactive_1d( &
      species, config, solution, changed, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    type(amr_multipatch_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: changed, ok

    type(amr_multipatch_reactive_solution_1d) :: backup
    type(amr_patch_set_1d) :: new_hierarchy
    type(amr_level_field_1d), allocatable :: old_fields(:), new_fields(:)
    type(amr_regrid_plan_collection_1d) :: collection
    logical :: local_ok
    integer :: transferred_cells, nx

    changed = .false.
    ok = .false.
    if (.not. solution%is_valid() .or. &
        .not. valid_tagged_multipatch_configuration( &
          config, size(solution%coarse, 1))) return
    backup = solution
    nx = solution%hierarchy%coarse_cells
    call build_tagged_patch_collection( &
      config, solution%coarse(:, 1:nx), collection, local_ok)
    if (.not. local_ok) return
    changed = .not. same_patch_plan(solution%hierarchy, collection)
    if (.not. changed) then
      solution%regrid_evaluations = backup%regrid_evaluations + 1
      ok = .true.
      return
    end if

    call extract_patch_fields(solution, old_fields, local_ok)
    if (.not. local_ok) return
    transferred_cells = overlapping_fine_cell_count( &
      solution%hierarchy, collection, config%amr_refinement_ratio)
    call regrid_patch_set_state_1d( &
      solution%coarse(:, 1:nx), solution%hierarchy, old_fields, collection, &
      config%amr_refinement_ratio, config%x_lower, config%x_upper, &
      new_hierarchy, new_fields, local_ok)
    if (.not. local_ok) then
      solution = backup
      return
    end if
    call install_patch_fields( &
      species, new_hierarchy, new_fields, solution, local_ok)
    if (.not. local_ok) then
      solution = backup
      return
    end if
    call recover_level_temperatures_1d( &
      species, solution%coarse, solution%coarse_temperature, nx, local_ok)
    if (.not. local_ok) then
      solution = backup
      return
    end if
    call refresh_multipatch_ghosts(species, config, solution, local_ok)
    if (.not. local_ok) then
      solution = backup
      return
    end if
    solution%regrid_evaluations = backup%regrid_evaluations + 1
    solution%regrids = backup%regrids + 1
    solution%overlap_cells_transferred = &
      backup%overlap_cells_transferred + transferred_cells
    ok = solution%is_valid()
    if (.not. ok) solution = backup
  end subroutine regrid_multipatch_reactive_1d

  subroutine simulate_multipatch_reactive_1d( &
      species, reactions, config, solution, initial_integrals, &
      final_integrals, ok, transport)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(reactive_1d_config), intent(in) :: config
    type(amr_multipatch_reactive_solution_1d), intent(out) :: solution
    real(dp), intent(out) :: initial_integrals(5), final_integrals(5)
    logical, intent(out) :: ok
    type(gas_transport_species), intent(in), optional :: transport(:)

    real(dp), allocatable :: all_integrals(:)
    real(dp) :: dt, tolerance
    logical :: local_ok, changed

    initial_integrals = 0.0_dp
    final_integrals = 0.0_dp
    ok = .false.
    if (config%transport_enabled .and. .not. present(transport)) return
    call initialize_tagged_multipatch_reactive_1d( &
      species, config, solution, local_ok)
    if (.not. local_ok) return
    allocate(all_integrals(reactive_nvar(size(species))))
    call multipatch_reactive_integrals_1d(solution, all_integrals, local_ok)
    if (.not. local_ok) return
    initial_integrals = all_integrals([irho, imx, imy, imz, iet])
    tolerance = 50.0_dp * epsilon(1.0_dp) * &
      max(1.0_dp, config%final_time)
    do while (solution%time < config%final_time - tolerance)
      if (solution%steps >= config%maximum_steps) return
      if (config%transport_enabled) then
        call multipatch_reactive_timestep_1d( &
          species, config, solution, dt, local_ok, transport)
      else
        call multipatch_reactive_timestep_1d( &
          species, config, solution, dt, local_ok)
      end if
      if (.not. local_ok) return
      dt = min(dt, config%final_time - solution%time)
      if (config%transport_enabled) then
        call advance_multipatch_reactive_1d( &
          species, reactions, config, dt, solution, local_ok, transport)
      else
        call advance_multipatch_reactive_1d( &
          species, reactions, config, dt, solution, local_ok)
      end if
      if (.not. local_ok) return
      if (mod(solution%steps, config%amr_regrid_interval) == 0) then
        call regrid_multipatch_reactive_1d( &
          species, config, solution, changed, local_ok)
        if (.not. local_ok) return
      end if
    end do
    solution%time = config%final_time
    call multipatch_reactive_integrals_1d(solution, all_integrals, local_ok)
    if (.not. local_ok) return
    final_integrals = all_integrals([irho, imx, imy, imz, iet])
    ok = .true.
  end subroutine simulate_multipatch_reactive_1d

  subroutine multipatch_reactive_integrals_1d(solution, integrals, ok)
    type(amr_multipatch_reactive_solution_1d), intent(in) :: solution
    real(dp), intent(out) :: integrals(:)
    logical, intent(out) :: ok

    type(amr_level_field_1d), allocatable :: fields(:)

    ok = solution%is_valid()
    if (.not. ok) return
    call extract_patch_fields(solution, fields, ok)
    if (.not. ok) return
    call composite_integral_patch_set_1d( &
      solution%coarse(:, 1:solution%hierarchy%coarse_cells), &
      fields, solution%hierarchy, integrals, ok)
  end subroutine multipatch_reactive_integrals_1d

  subroutine write_multipatch_reactive_1d_csv(path, species, solution, ok)
    character(len=*), intent(in) :: path
    type(nasa7_species), intent(in) :: species(:)
    type(amr_multipatch_reactive_solution_1d), intent(in) :: solution
    logical, intent(out) :: ok

    real(dp), allocatable :: q(:)
    real(dp) :: x, patch_lower
    logical :: local_ok
    integer :: unit, status, k, patch, cell, next_parent, fine_cells
    integer :: lower, upper

    ok = .false.
    if (.not. solution%is_valid()) return
    allocate(q(reactive_nprim(size(species))))
    open(newunit=unit, file=trim(path), status="replace", action="write", &
      iostat=status)
    if (status /= 0) return
    write(unit, '(a)', advance='no') &
      "level,cell_dx,time,x,rho,u,v,w,pressure,temperature,rhoE"
    do k = 1, size(species)
      write(unit, '(a)', advance='no') ",Y_" // trim(species(k)%name)
    end do
    write(unit, '(a)') ""

    next_parent = 1
    do patch = 1, solution%patch_count()
      lower = solution%hierarchy%patches(patch)%fine_coarse_lower
      upper = solution%hierarchy%patches(patch)%fine_coarse_upper
      do cell = next_parent, lower - 1
        x = solution%hierarchy%x_lower + &
          (real(cell, dp) - 0.5_dp) * solution%hierarchy%coarse_dx
        call write_amr_cell( &
          unit, 0, solution%hierarchy%coarse_dx, solution%time, x, species, &
          solution%coarse(:, cell), solution%coarse_temperature(cell), &
          q, local_ok)
        if (.not. local_ok) then
          close(unit)
          return
        end if
      end do
      patch_lower = solution%hierarchy%x_lower + &
        real(lower - 1, dp) * solution%hierarchy%coarse_dx
      fine_cells = solution%hierarchy%patches(patch)%fine%cell_count()
      do cell = 1, fine_cells
        x = patch_lower + &
          (real(cell, dp) - 0.5_dp) * solution%hierarchy%fine_dx
        call write_amr_cell( &
          unit, 1, solution%hierarchy%fine_dx, solution%time, x, species, &
          solution%patches(patch)%state(:, cell), &
          solution%patches(patch)%temperature(cell), q, local_ok)
        if (.not. local_ok) then
          close(unit)
          return
        end if
      end do
      next_parent = upper + 1
    end do
    do cell = next_parent, solution%hierarchy%coarse_cells
      x = solution%hierarchy%x_lower + &
        (real(cell, dp) - 0.5_dp) * solution%hierarchy%coarse_dx
      call write_amr_cell( &
        unit, 0, solution%hierarchy%coarse_dx, solution%time, x, species, &
        solution%coarse(:, cell), solution%coarse_temperature(cell), &
        q, local_ok)
      if (.not. local_ok) then
        close(unit)
        return
      end if
    end do
    close(unit)
    ok = .true.
  end subroutine write_multipatch_reactive_1d_csv

  subroutine extract_patch_fields(solution, fields, ok)
    type(amr_multipatch_reactive_solution_1d), intent(in) :: solution
    type(amr_level_field_1d), allocatable, intent(out) :: fields(:)
    logical, intent(out) :: ok

    integer :: patch, fine_cells, nvar

    ok = solution%is_valid()
    if (.not. ok) return
    nvar = size(solution%coarse, 1)
    allocate(fields(solution%patch_count()))
    do patch = 1, size(fields)
      fine_cells = solution%hierarchy%patches(patch)%fine%cell_count()
      allocate(fields(patch)%values(nvar, fine_cells))
      fields(patch)%values = &
        solution%patches(patch)%state(:, 1:fine_cells)
    end do
    ok = .true.
  end subroutine extract_patch_fields

  subroutine install_patch_fields( &
      species, hierarchy, fields, solution, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_set_1d), intent(in) :: hierarchy
    type(amr_level_field_1d), intent(in) :: fields(:)
    type(amr_multipatch_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok

    logical :: local_ok
    integer :: patch, fine_cells, nvar

    ok = hierarchy%is_valid() .and. size(fields) == hierarchy%patch_count()
    if (.not. ok) return
    nvar = size(solution%coarse, 1)
    solution%hierarchy = hierarchy
    if (allocated(solution%patches)) deallocate(solution%patches)
    allocate(solution%patches(hierarchy%patch_count()))
    do patch = 1, hierarchy%patch_count()
      fine_cells = hierarchy%patches(patch)%fine%cell_count()
      if (.not. allocated(fields(patch)%values) .or. &
          size(fields(patch)%values, 1) /= nvar .or. &
          size(fields(patch)%values, 2) /= fine_cells) then
        ok = .false.
        return
      end if
      allocate(solution%patches(patch)%state(nvar, 0:fine_cells + 1))
      allocate(solution%patches(patch)%temperature(0:fine_cells + 1))
      allocate(solution%patches(patch)%left_ghost_state( &
        nvar, amr_ppm_ghost_width))
      allocate(solution%patches(patch)%right_ghost_state( &
        nvar, amr_ppm_ghost_width))
      allocate(solution%patches(patch)%left_ghost_temperature( &
        amr_ppm_ghost_width))
      allocate(solution%patches(patch)%right_ghost_temperature( &
        amr_ppm_ghost_width))
      solution%patches(patch)%state = 0.0_dp
      solution%patches(patch)%temperature = 0.0_dp
      solution%patches(patch)%left_ghost_state = 0.0_dp
      solution%patches(patch)%right_ghost_state = 0.0_dp
      solution%patches(patch)%left_ghost_temperature = 0.0_dp
      solution%patches(patch)%right_ghost_temperature = 0.0_dp
      solution%patches(patch)%state(:, 1:fine_cells) = fields(patch)%values
      call recover_level_temperatures_1d( &
        species, solution%patches(patch)%state, &
        solution%patches(patch)%temperature, fine_cells, local_ok)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
    end do
    ok = .true.
  end subroutine install_patch_fields

  subroutine build_tagged_patch_collection( &
      config, state, collection, ok)
    type(reactive_1d_config), intent(in) :: config
    real(dp), intent(in) :: state(:, :)
    type(amr_regrid_plan_collection_1d), intent(out) :: collection
    logical, intent(out) :: ok

    type(amr_tagging_criteria_1d) :: criteria
    logical, allocatable :: tags(:)
    integer :: nx, patch, parent_buffer, allowed_lower, allowed_upper

    nx = size(state, 2)
    ok = valid_tagged_multipatch_configuration(config, size(state, 1)) .and. &
      nx == config%nx
    if (.not. ok) return
    call multipatch_criteria_from_config(config, criteria)
    allocate(tags(nx))
    call tag_gradient_1d(state, criteria, tags, ok)
    if (.not. ok) return

    parent_buffer = 1
    if (uses_ppm_reconstruction(config)) then
      parent_buffer = (amr_ppm_ghost_width + &
        config%amr_refinement_ratio - 1) / &
        config%amr_refinement_ratio + 1
    end if
    if (trim(config%boundary_condition) == "outflow") then
      allowed_lower = 1
      allowed_upper = nx
    else
      allowed_lower = parent_buffer + 1
      allowed_upper = nx - parent_buffer
    end if
    if (allowed_lower > allowed_upper .or. &
        allowed_upper - allowed_lower + 1 < criteria%minimum_patch_cells) then
      tags = .false.
    else
      if (allowed_lower > 1) tags(1:allowed_lower - 1) = .false.
      if (allowed_upper < nx) tags(allowed_upper + 1:nx) = .false.
    end if
    call build_regrid_plan_collection_1d( &
      tags, criteria%buffer_cells, criteria%minimum_patch_cells, &
      criteria%maximum_patch_gap_cells, collection, ok)
    if (.not. ok) return
    do patch = 1, collection%patch_count()
      collection%plans(patch)%patch_lower = max( &
        collection%plans(patch)%patch_lower, allowed_lower)
      collection%plans(patch)%patch_upper = min( &
        collection%plans(patch)%patch_upper, allowed_upper)
      do while (collection%plans(patch)%patch_upper - &
          collection%plans(patch)%patch_lower + 1 < &
          criteria%minimum_patch_cells)
        if (collection%plans(patch)%patch_lower > allowed_lower) &
          collection%plans(patch)%patch_lower = &
            collection%plans(patch)%patch_lower - 1
        if (collection%plans(patch)%patch_upper - &
            collection%plans(patch)%patch_lower + 1 >= &
            criteria%minimum_patch_cells) exit
        if (collection%plans(patch)%patch_upper < allowed_upper) &
          collection%plans(patch)%patch_upper = &
            collection%plans(patch)%patch_upper + 1
      end do
    end do
    ok = collection%is_valid()
  end subroutine build_tagged_patch_collection

  pure logical function same_patch_plan(hierarchy, collection) result(same)
    type(amr_patch_set_1d), intent(in) :: hierarchy
    type(amr_regrid_plan_collection_1d), intent(in) :: collection

    integer :: patch

    same = hierarchy%is_valid() .and. collection%is_valid() .and. &
      hierarchy%coarse_cells == collection%coarse_cells .and. &
      hierarchy%patch_count() == collection%patch_count()
    if (.not. same) return
    do patch = 1, hierarchy%patch_count()
      same = hierarchy%patches(patch)%fine_coarse_lower == &
        collection%plans(patch)%patch_lower .and. &
        hierarchy%patches(patch)%fine_coarse_upper == &
        collection%plans(patch)%patch_upper
      if (.not. same) return
    end do
  end function same_patch_plan

  pure integer function overlapping_fine_cell_count( &
      old_hierarchy, collection, refinement_ratio) result(count)
    type(amr_patch_set_1d), intent(in) :: old_hierarchy
    type(amr_regrid_plan_collection_1d), intent(in) :: collection
    integer, intent(in) :: refinement_ratio

    integer :: old_patch, new_patch, overlap_lower, overlap_upper

    count = 0
    if (.not. old_hierarchy%is_valid() .or. &
        .not. collection%is_valid() .or. refinement_ratio < 1) return
    if (old_hierarchy%refinement_ratio /= refinement_ratio) return
    do new_patch = 1, collection%patch_count()
      do old_patch = 1, old_hierarchy%patch_count()
        overlap_lower = max( &
          collection%plans(new_patch)%patch_lower, &
          old_hierarchy%patches(old_patch)%fine_coarse_lower)
        overlap_upper = min( &
          collection%plans(new_patch)%patch_upper, &
          old_hierarchy%patches(old_patch)%fine_coarse_upper)
        if (overlap_lower <= overlap_upper) count = count + &
          (overlap_upper - overlap_lower + 1) * refinement_ratio
      end do
    end do
  end function overlapping_fine_cell_count

  pure subroutine multipatch_criteria_from_config(config, criteria)
    type(reactive_1d_config), intent(in) :: config
    type(amr_tagging_criteria_1d), intent(out) :: criteria

    criteria%component = config%amr_tag_component
    criteria%relative_gradient_threshold = &
      config%amr_relative_gradient_threshold
    criteria%absolute_gradient_threshold = &
      config%amr_absolute_gradient_threshold
    criteria%scale_floor = config%amr_scale_floor
    criteria%buffer_cells = config%amr_buffer_cells
    criteria%minimum_patch_cells = config%amr_minimum_patch_cells
    criteria%maximum_patch_gap_cells = config%amr_maximum_patch_gap_cells
  end subroutine multipatch_criteria_from_config

  pure logical function valid_tagged_multipatch_configuration( &
      config, nvar) result(valid)
    type(reactive_1d_config), intent(in) :: config
    integer, intent(in) :: nvar

    valid = config%amr_enabled .and. config%amr_multipatch_enabled .and. &
      config%amr_max_levels == 2 .and. config%nx >= 8 .and. &
      config%amr_refinement_ratio >= 2 .and. &
      config%amr_regrid_interval >= 1 .and. &
      config%amr_tag_component >= 1 .and. &
      config%amr_tag_component <= nvar .and. &
      config%amr_buffer_cells >= 0 .and. &
      config%amr_minimum_patch_cells >= 1 .and. &
      config%amr_minimum_patch_cells <= config%nx - 2 .and. &
      config%amr_maximum_patch_gap_cells >= 0 .and. &
      config%amr_relative_gradient_threshold >= 0.0_dp .and. &
      config%amr_absolute_gradient_threshold >= 0.0_dp .and. &
      config%amr_scale_floor > 0.0_dp .and. &
      (trim(config%amr_reconstruction) == "pcm" .or. &
        trim(config%amr_reconstruction) == "plm" .or. &
        trim(config%amr_reconstruction) == "ppm" .or. &
        trim(config%amr_reconstruction) == "characteristic_ppm")
  end function valid_tagged_multipatch_configuration

  subroutine refresh_multipatch_ghosts(species, config, solution, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_1d_config), intent(in) :: config
    type(amr_multipatch_reactive_solution_1d), intent(inout) :: solution
    logical, intent(out) :: ok

    logical :: local_ok
    integer :: patch

    call fill_physical_ghosts_1d( &
      solution%coarse, solution%coarse_temperature, &
      solution%hierarchy%coarse_cells, config%boundary_condition, local_ok)
    if (.not. local_ok) then
      ok = .false.
      return
    end if
    do patch = 1, solution%patch_count()
      call fill_fine_ghosts_1d( &
        species, solution%hierarchy%patches(patch), &
        solution%coarse, solution%coarse, 1.0_dp, &
        solution%patches(patch)%state, &
        solution%patches(patch)%temperature, local_ok, &
        config%boundary_condition)
      if (.not. local_ok) then
        ok = .false.
        return
      end if
      if (uses_ppm_reconstruction(config)) then
        call fill_fine_wide_ghosts_1d( &
          species, solution%hierarchy%patches(patch), &
          solution%coarse, solution%coarse, 1.0_dp, &
          solution%patches(patch)%left_ghost_state, &
          solution%patches(patch)%right_ghost_state, &
          solution%patches(patch)%left_ghost_temperature, &
          solution%patches(patch)%right_ghost_temperature, local_ok, &
          solution%patches(patch)%state, &
          solution%patches(patch)%temperature, config%boundary_condition)
        if (.not. local_ok) then
          ok = .false.
          return
        end if
      end if
    end do
    ok = .true.
  end subroutine refresh_multipatch_ghosts

  pure logical function uses_ppm_reconstruction(config) result(enabled)
    type(reactive_1d_config), intent(in) :: config

    enabled = trim(config%amr_reconstruction) == "ppm" .or. &
      trim(config%amr_reconstruction) == "characteristic_ppm"
  end function uses_ppm_reconstruction

end module amr_multipatch_reactive_1d_mod
