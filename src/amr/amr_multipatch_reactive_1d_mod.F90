module amr_multipatch_reactive_1d_mod
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use transport_database_mod, only: gas_transport_species
  use simulation_config_reactive_1d_mod, only: reactive_1d_config
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_cfl_timestep, reactive_transport_timestep, &
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
    fill_fine_ghosts_1d, fill_fine_wide_ghosts_1d
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
  contains
    procedure :: patch_count => multipatch_reactive_patch_count
    procedure :: is_valid => multipatch_reactive_is_valid
  end type amr_multipatch_reactive_solution_1d

  public :: initialize_multipatch_reactive_1d
  public :: multipatch_reactive_timestep_1d
  public :: advance_multipatch_reactive_1d
  public :: advance_multipatch_reactive_hydro_1d
  public :: multipatch_reactive_integrals_1d

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
      self%steps >= 0
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
    call refresh_multipatch_ghosts(species, config, solution, local_ok)
    if (.not. local_ok) return
    ok = solution%is_valid() .and. coarse_dx > 0.0_dp
  end subroutine initialize_multipatch_reactive_1d

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
