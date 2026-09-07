module amr_reactive_3d_checkpoint_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: &
    elementary_reaction, valid_elementary_reaction, &
    valid_arrhenius_rate, valid_troe_parameters
  use gas_transport_mod, only: &
    gas_transport_species, compatible_transport_database
  use reactive_1d_mod, only: reactive_nvar
  use reactive_3d_mod, only: recover_reactive_temperatures_3d
  use simulation_config_reactive_3d_mod, only: reactive_3d_config
  use simulation_config_amr_reactive_3d_mod, only: amr_reactive_3d_config
  use amr_hierarchy_3d_mod, only: &
    amr_patch_3d, restrict_average_3d
  implicit none
  private

  character(len=*), parameter :: checkpoint_magic = &
    "PELEF_AMR_REACTIVE_3D_CHECKPOINT"
  integer, parameter :: checkpoint_schema_fixed = 2
  integer, parameter :: checkpoint_schema_selected = 3
  integer, parameter :: checkpoint_schema_fixed_transport = 4
  integer, parameter :: checkpoint_schema_selected_transport = 5
  integer, parameter :: bundle_sha256_length = 64
  character(len=*), parameter, public :: &
    amr_reactive_3d_transport_operator = "STATIC_AMR_3D_R_T_H_T_R_V1"
  character(len=*), parameter :: transport_parameter_convention = &
    "EPSILON_OVER_K_K;SIGMA_ANGSTROM;DIPOLE_DEBYE;POLARIZABILITY_ANGSTROM3;ROT_RELAX_DIMENSIONLESS"
  character(len=*), parameter :: transport_checkpoint_phase = &
    "POST_ACCEPTED_COARSE_STEP"

  public :: write_amr_reactive_3d_checkpoint
  public :: read_amr_reactive_3d_checkpoint

contains

  subroutine write_amr_reactive_3d_checkpoint( &
      path, species, flow_config, amr_config, patch, &
      coarse_state, coarse_temperature, fine_state, fine_temperature, &
      time, steps, initial_integrals, maximum_reflux, ok, message, &
      bundle_sha256, chemistry_integrator, base_mole_fractions, reactions, &
      transport, maximum_transport_diffusivity, minimum_transport_theta)
    character(len=*), intent(in) :: path
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_3d_config), intent(in) :: flow_config
    type(amr_reactive_3d_config), intent(in) :: amr_config
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: coarse_state(:, :, :, :)
    real(dp), intent(in) :: coarse_temperature(:, :, :)
    real(dp), intent(in) :: fine_state(:, :, :, :)
    real(dp), intent(in) :: fine_temperature(:, :, :)
    real(dp), intent(in) :: time
    integer, intent(in) :: steps
    real(dp), intent(in) :: initial_integrals(:), maximum_reflux
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    character(len=*), intent(in), optional :: bundle_sha256
    character(len=*), intent(in), optional :: chemistry_integrator
    real(dp), intent(in), optional :: base_mole_fractions(:)
    type(elementary_reaction), intent(in), optional :: reactions(:)
    type(gas_transport_species), intent(in), optional :: transport(:)
    real(dp), intent(in), optional :: maximum_transport_diffusivity
    real(dp), intent(in), optional :: minimum_transport_theta

    integer :: unit, status, i, j, k, species_index, nvar, schema
    logical :: synchronized, temperature_consistent, selected_context
    logical :: transport_context

    ok = .false.
    message = ""
    nvar = reactive_nvar(size(species))
    selected_context = present(bundle_sha256) .or. &
      present(chemistry_integrator) .or. present(base_mole_fractions) .or. &
      present(reactions)
    transport_context = present(transport) .or. &
      present(maximum_transport_diffusivity) .or. &
      present(minimum_transport_theta)
    if (len_trim(path) == 0) then
      message = "3D AMR checkpoint path is empty"
      return
    end if
    if (selected_context) then
      if (.not. present(bundle_sha256) .or. &
          .not. present(chemistry_integrator) .or. &
          .not. present(base_mole_fractions) .or. &
          .not. present(reactions)) then
        message = "Selected 3D AMR checkpoint context is incomplete"
        return
      end if
      if (.not. valid_selected_checkpoint_context( &
          size(species), bundle_sha256, chemistry_integrator, &
          base_mole_fractions, reactions, &
          flow_config%chemistry_relative_tolerance, &
          flow_config%chemistry_absolute_tolerance)) then
        message = "Selected 3D AMR checkpoint context is invalid"
        return
      end if
    end if
    if (flow_config%transport_enabled) then
      if (.not. present(transport) .or. &
          .not. present(maximum_transport_diffusivity) .or. &
          .not. present(minimum_transport_theta)) then
        message = "Transport 3D AMR checkpoint context is incomplete"
        return
      end if
      if (.not. valid_transport_checkpoint_context( &
          species, flow_config, transport, maximum_transport_diffusivity, &
          minimum_transport_theta)) then
        message = "Transport 3D AMR checkpoint context is invalid"
        return
      end if
    else if (transport_context) then
      message = "Transport 3D AMR checkpoint context is unexpected"
      return
    end if
    if (selected_context) then
      if (flow_config%transport_enabled) then
        schema = checkpoint_schema_selected_transport
      else
        schema = checkpoint_schema_selected
      end if
    else
      if (flow_config%transport_enabled) then
        schema = checkpoint_schema_fixed_transport
      else
        schema = checkpoint_schema_fixed
      end if
    end if
    if (.not. valid_checkpoint_shapes( &
        species, patch, coarse_state, coarse_temperature, &
        fine_state, fine_temperature)) then
      message = "Invalid 3D AMR checkpoint state"
      return
    end if
    if (size(initial_integrals) /= nvar) then
      message = "Invalid 3D AMR checkpoint state"
      return
    end if
    if (.not. all(ieee_is_finite(initial_integrals)) .or. &
        .not. all(ieee_is_finite(coarse_state)) .or. &
        .not. all(ieee_is_finite(coarse_temperature)) .or. &
        .not. all(ieee_is_finite(fine_state)) .or. &
        .not. all(ieee_is_finite(fine_temperature))) then
      message = "Invalid 3D AMR checkpoint state"
      return
    end if
    if (minval(coarse_temperature) <= 0.0_dp .or. &
        minval(fine_temperature) <= 0.0_dp) then
      message = "Invalid 3D AMR checkpoint state"
      return
    end if
    if (.not. ieee_is_finite(time) .or. &
        .not. ieee_is_finite(maximum_reflux) .or. &
        .not. ieee_is_finite(flow_config%final_time)) then
      message = "Invalid 3D AMR checkpoint state"
      return
    end if
    if (time < 0.0_dp .or. maximum_reflux < 0.0_dp .or. steps < 0 .or. &
        time > flow_config%final_time + &
          64.0_dp * epsilon(1.0_dp) * &
            max(1.0_dp, abs(flow_config%final_time)) .or. &
        steps > flow_config%maximum_steps .or. &
        (steps == 0 .and. time > tiny(1.0_dp)) .or. &
        (steps > 0 .and. time <= 0.0_dp)) then
      message = "Invalid 3D AMR checkpoint state"
      return
    end if
    call checkpoint_levels_synchronized( &
      patch, coarse_state, fine_state, synchronized)
    if (.not. synchronized) then
      message = "3D AMR checkpoint levels are not synchronized"
      return
    end if
    call checkpoint_temperature_is_consistent( &
      species, coarse_state, coarse_temperature, temperature_consistent)
    if (.not. temperature_consistent) then
      message = "3D AMR coarse checkpoint temperature is inconsistent"
      return
    end if
    call checkpoint_temperature_is_consistent( &
      species, fine_state, fine_temperature, temperature_consistent)
    if (.not. temperature_consistent) then
      message = "3D AMR fine checkpoint temperature is inconsistent"
      return
    end if
    if (.not. checkpoint_configuration_matches( &
          flow_config, amr_config, patch, selected_context)) then
      message = "Invalid 3D AMR checkpoint configuration"
      return
    end if

    open(newunit=unit, file=trim(path), status="replace", action="write", &
      form="formatted", iostat=status)
    if (status /= 0) then
      write(message, '(a,1x,a)') &
        "Could not create 3D AMR checkpoint:", trim(path)
      return
    end if
    write(unit, '(a)', iostat=status) checkpoint_magic
    if (status /= 0) go to 900
    write(unit, '(*(i0,1x))', iostat=status) &
      schema, size(species), nvar
    if (status /= 0) go to 900
    if (selected_context) then
      call write_selected_checkpoint_context( &
        unit, reactions, bundle_sha256, chemistry_integrator, &
        base_mole_fractions, flow_config%chemistry_relative_tolerance, &
        flow_config%chemistry_absolute_tolerance, status)
      if (status /= 0) go to 900
    end if
    if (flow_config%transport_enabled) then
      call write_transport_checkpoint_context( &
        unit, transport, flow_config, maximum_transport_diffusivity, &
        minimum_transport_theta, status)
      if (status /= 0) go to 900
    end if
    do species_index = 1, size(species)
      write(unit, '(a)', iostat=status) trim(species(species_index)%name)
      if (status /= 0) go to 900
      write(unit, '(*(es27.18e3,1x))', iostat=status) &
        species(species_index)%molecular_weight, &
        species(species_index)%temperature_min, &
        species(species_index)%temperature_mid, &
        species(species_index)%temperature_max, &
        species(species_index)%low_coefficients, &
        species(species_index)%high_coefficients
      if (status /= 0) go to 900
    end do
    write(unit, '(a)', iostat=status) trim(flow_config%thermo_model)
    if (status /= 0) go to 900
    write(unit, '(a)', iostat=status) trim(flow_config%riemann_solver)
    if (status /= 0) go to 900
    write(unit, '(a)', iostat=status) trim(flow_config%boundary_condition)
    if (status /= 0) go to 900
    write(unit, '(a)', iostat=status) trim(flow_config%reconstruction)
    if (status /= 0) go to 900
    write(unit, '(a)', iostat=status) trim(flow_config%limiter)
    if (status /= 0) go to 900
    write(unit, '(*(i0,1x))', iostat=status) &
      flow_config%nx, flow_config%ny, flow_config%nz, &
      patch%coarse_i_lower, patch%coarse_i_upper, &
      patch%coarse_j_lower, patch%coarse_j_upper, &
      patch%coarse_k_lower, patch%coarse_k_upper, patch%refinement_ratio
    if (status /= 0) go to 900
    write(unit, '(*(es27.18e3,1x))', iostat=status) &
      flow_config%x_lower, flow_config%x_upper, &
      flow_config%y_lower, flow_config%y_upper, &
      flow_config%z_lower, flow_config%z_upper, flow_config%cfl
    if (status /= 0) go to 900
    write(unit, '(*(i0,1x))', iostat=status) &
      merge(1, 0, flow_config%chemistry_enabled), &
      merge(1, 0, flow_config%transport_enabled)
    if (status /= 0) go to 900
    write(unit, '(2(es27.18e3,1x),i0)', iostat=status) &
      time, maximum_reflux, steps
    if (status /= 0) go to 900
    write(unit, '(*(es27.18e3,1x))', iostat=status) initial_integrals
    if (status /= 0) go to 900
    write(unit, '(*(i0,1x))', iostat=status) &
      size(coarse_state, 2), size(coarse_state, 3), size(coarse_state, 4)
    if (status /= 0) go to 900
    do k = 1, size(coarse_state, 4)
      do j = 1, size(coarse_state, 3)
        do i = 1, size(coarse_state, 2)
          write(unit, '(*(es27.18e3,1x))', iostat=status) &
            coarse_state(:, i, j, k), coarse_temperature(i, j, k)
          if (status /= 0) go to 900
        end do
      end do
    end do
    write(unit, '(*(i0,1x))', iostat=status) &
      size(fine_state, 2), size(fine_state, 3), size(fine_state, 4)
    if (status /= 0) go to 900
    do k = 1, size(fine_state, 4)
      do j = 1, size(fine_state, 3)
        do i = 1, size(fine_state, 2)
          write(unit, '(*(es27.18e3,1x))', iostat=status) &
            fine_state(:, i, j, k), fine_temperature(i, j, k)
          if (status /= 0) go to 900
        end do
      end do
    end do
    write(unit, '(a)', iostat=status) "END_CHECKPOINT"
    if (status /= 0) go to 900
    close(unit, iostat=status)
    if (status /= 0) then
      message = "Failed to close 3D AMR checkpoint"
      return
    end if
    ok = .true.
    return

900 continue
    close(unit)
    message = "Failed while writing 3D AMR checkpoint"
  end subroutine write_amr_reactive_3d_checkpoint

  subroutine read_amr_reactive_3d_checkpoint( &
      path, species, flow_config, amr_config, patch, &
      coarse_state, coarse_temperature, fine_state, fine_temperature, &
      time, steps, initial_integrals, maximum_reflux, ok, message, &
      bundle_sha256, chemistry_integrator, base_mole_fractions, reactions, &
      transport, maximum_transport_diffusivity, minimum_transport_theta)
    character(len=*), intent(in) :: path
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_3d_config), intent(in) :: flow_config
    type(amr_reactive_3d_config), intent(in) :: amr_config
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(inout) :: coarse_state(:, :, :, :)
    real(dp), intent(inout) :: coarse_temperature(:, :, :)
    real(dp), intent(inout) :: fine_state(:, :, :, :)
    real(dp), intent(inout) :: fine_temperature(:, :, :)
    real(dp), intent(inout) :: time
    integer, intent(inout) :: steps
    real(dp), intent(inout) :: initial_integrals(:), maximum_reflux
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    character(len=*), intent(in), optional :: bundle_sha256
    character(len=*), intent(in), optional :: chemistry_integrator
    real(dp), intent(in), optional :: base_mole_fractions(:)
    type(elementary_reaction), intent(in), optional :: reactions(:)
    type(gas_transport_species), intent(in), optional :: transport(:)
    real(dp), intent(inout), optional :: maximum_transport_diffusivity
    real(dp), intent(inout), optional :: minimum_transport_theta

    real(dp), allocatable :: candidate_coarse_state(:, :, :, :)
    real(dp), allocatable :: candidate_coarse_temperature(:, :, :)
    real(dp), allocatable :: candidate_fine_state(:, :, :, :)
    real(dp), allocatable :: candidate_fine_temperature(:, :, :)
    real(dp), allocatable :: recovered_coarse_temperature(:, :, :)
    real(dp), allocatable :: recovered_fine_temperature(:, :, :)
    real(dp), allocatable :: candidate_integrals(:)
    real(dp) :: stored_species_values(18), expected_species_values(18)
    real(dp) :: stored_mesh(7), expected_mesh(7)
    real(dp) :: stored_time, stored_maximum_reflux, time_tolerance
    real(dp) :: candidate_maximum_transport_diffusivity
    real(dp) :: candidate_minimum_transport_theta
    character(len=1024) :: magic, stored_name, stored_thermo
    character(len=1024) :: stored_solver, stored_boundary, end_marker
    character(len=1024) :: stored_reconstruction, stored_limiter
    integer :: stored_header(3), stored_layout(10), expected_layout(10)
    integer :: stored_flags(2), expected_flags(2), stored_coarse_shape(3)
    integer :: stored_fine_shape(3), stored_steps
    integer :: unit, status, i, j, k, species_index, nvar, expected_schema
    logical :: local_ok, synchronized, temperature_consistent
    logical :: selected_context, context_ok, transport_context

    ok = .false.
    message = ""
    candidate_maximum_transport_diffusivity = 0.0_dp
    candidate_minimum_transport_theta = 1.0_dp
    nvar = reactive_nvar(size(species))
    selected_context = present(bundle_sha256) .or. &
      present(chemistry_integrator) .or. present(base_mole_fractions) .or. &
      present(reactions)
    transport_context = present(transport) .or. &
      present(maximum_transport_diffusivity) .or. &
      present(minimum_transport_theta)
    if (len_trim(path) == 0) then
      message = "3D AMR restart path is empty"
      return
    end if
    if (selected_context) then
      if (.not. present(bundle_sha256) .or. &
          .not. present(chemistry_integrator) .or. &
          .not. present(base_mole_fractions) .or. &
          .not. present(reactions)) then
        message = "Selected 3D AMR restart context is incomplete"
        return
      end if
      if (.not. valid_selected_checkpoint_context( &
          size(species), bundle_sha256, chemistry_integrator, &
          base_mole_fractions, reactions, &
          flow_config%chemistry_relative_tolerance, &
          flow_config%chemistry_absolute_tolerance)) then
        message = "Selected 3D AMR restart context is invalid"
        return
      end if
    end if
    if (flow_config%transport_enabled) then
      if (.not. present(transport) .or. &
          .not. present(maximum_transport_diffusivity) .or. &
          .not. present(minimum_transport_theta)) then
        message = "Transport 3D AMR restart context is incomplete"
        return
      end if
      if (.not. compatible_transport_database(species, transport)) then
        message = "Transport 3D AMR restart context is invalid"
        return
      end if
    else if (transport_context) then
      message = "Transport 3D AMR restart context is unexpected"
      return
    end if
    if (selected_context) then
      if (flow_config%transport_enabled) then
        expected_schema = checkpoint_schema_selected_transport
      else
        expected_schema = checkpoint_schema_selected
      end if
    else
      if (flow_config%transport_enabled) then
        expected_schema = checkpoint_schema_fixed_transport
      else
        expected_schema = checkpoint_schema_fixed
      end if
    end if
    if (.not. valid_checkpoint_shapes( &
          species, patch, coarse_state, coarse_temperature, &
          fine_state, fine_temperature) .or. &
        size(initial_integrals) /= nvar .or. &
        .not. checkpoint_configuration_matches( &
          flow_config, amr_config, patch, selected_context)) then
      message = "Invalid 3D AMR restart target"
      return
    end if
    open(newunit=unit, file=trim(path), status="old", action="read", &
      form="formatted", iostat=status)
    if (status /= 0) then
      write(message, '(a,1x,a)') &
        "Could not open 3D AMR checkpoint:", trim(path)
      return
    end if
    read(unit, '(a)', iostat=status) magic
    if (status /= 0) go to 900
    if (trim(magic) /= checkpoint_magic) go to 900
    read(unit, *, iostat=status) stored_header
    if (status /= 0) go to 900
    if (any(stored_header /= [expected_schema, size(species), nvar])) &
      go to 900
    if (selected_context) then
      call read_selected_checkpoint_context( &
        unit, reactions, bundle_sha256, chemistry_integrator, &
        base_mole_fractions, flow_config%chemistry_relative_tolerance, &
        flow_config%chemistry_absolute_tolerance, context_ok, message)
      if (.not. context_ok) go to 900
    end if
    if (flow_config%transport_enabled) then
      call read_transport_checkpoint_context( &
        unit, species, transport, flow_config, &
        candidate_maximum_transport_diffusivity, &
        candidate_minimum_transport_theta, context_ok, message)
      if (.not. context_ok) go to 900
    end if
    do species_index = 1, size(species)
      read(unit, '(a)', iostat=status) stored_name
      if (status /= 0) go to 900
      if (trim(stored_name) /= trim(species(species_index)%name)) go to 900
      read(unit, *, iostat=status) stored_species_values
      if (status /= 0) go to 900
      expected_species_values = [ &
        species(species_index)%molecular_weight, &
        species(species_index)%temperature_min, &
        species(species_index)%temperature_mid, &
        species(species_index)%temperature_max, &
        species(species_index)%low_coefficients, &
        species(species_index)%high_coefficients]
      if (.not. all(checkpoint_real_matches( &
          stored_species_values, expected_species_values))) go to 900
    end do
    read(unit, '(a)', iostat=status) stored_thermo
    if (status /= 0) go to 900
    if (trim(stored_thermo) /= trim(flow_config%thermo_model)) go to 900
    read(unit, '(a)', iostat=status) stored_solver
    if (status /= 0) go to 900
    if (trim(stored_solver) /= trim(flow_config%riemann_solver)) go to 900
    read(unit, '(a)', iostat=status) stored_boundary
    if (status /= 0) go to 900
    if (trim(stored_boundary) /= trim(flow_config%boundary_condition)) &
      go to 900
    read(unit, '(a)', iostat=status) stored_reconstruction
    if (status /= 0) go to 900
    if (trim(stored_reconstruction) /= trim(flow_config%reconstruction)) &
      go to 900
    read(unit, '(a)', iostat=status) stored_limiter
    if (status /= 0) go to 900
    if (trim(stored_limiter) /= trim(flow_config%limiter)) go to 900
    read(unit, *, iostat=status) stored_layout
    if (status /= 0) go to 900
    expected_layout = [ &
      flow_config%nx, flow_config%ny, flow_config%nz, &
      patch%coarse_i_lower, patch%coarse_i_upper, &
      patch%coarse_j_lower, patch%coarse_j_upper, &
      patch%coarse_k_lower, patch%coarse_k_upper, patch%refinement_ratio]
    if (any(stored_layout /= expected_layout)) go to 900
    read(unit, *, iostat=status) stored_mesh
    if (status /= 0) go to 900
    expected_mesh = [ &
      flow_config%x_lower, flow_config%x_upper, &
      flow_config%y_lower, flow_config%y_upper, &
      flow_config%z_lower, flow_config%z_upper, flow_config%cfl]
    if (.not. all(checkpoint_real_matches( &
        stored_mesh, expected_mesh))) go to 900
    read(unit, *, iostat=status) stored_flags
    if (status /= 0) go to 900
    expected_flags = [ &
      merge(1, 0, flow_config%chemistry_enabled), &
      merge(1, 0, flow_config%transport_enabled)]
    if (any(stored_flags /= expected_flags)) go to 900
    read(unit, *, iostat=status) &
      stored_time, stored_maximum_reflux, stored_steps
    if (status /= 0) go to 900
    if (.not. ieee_is_finite(stored_time) .or. &
        .not. ieee_is_finite(stored_maximum_reflux) .or. &
        .not. ieee_is_finite(flow_config%final_time)) go to 900
    time_tolerance = 64.0_dp * epsilon(1.0_dp) * &
      max(1.0_dp, abs(flow_config%final_time))
    if (stored_time < 0.0_dp .or. &
        stored_time > flow_config%final_time + time_tolerance .or. &
        stored_maximum_reflux < 0.0_dp .or. stored_steps < 0 .or. &
        stored_steps > flow_config%maximum_steps .or. &
        (stored_steps == 0 .and. stored_time > tiny(1.0_dp)) .or. &
        (stored_steps > 0 .and. stored_time <= 0.0_dp)) go to 900
    allocate(candidate_integrals(nvar))
    read(unit, *, iostat=status) candidate_integrals
    if (status /= 0) go to 900
    if (.not. all(ieee_is_finite(candidate_integrals))) go to 900

    read(unit, *, iostat=status) stored_coarse_shape
    if (status /= 0) go to 900
    if (any(stored_coarse_shape /= &
        [patch%coarse_nx, patch%coarse_ny, patch%coarse_nz])) go to 900
    allocate(candidate_coarse_state(nvar, patch%coarse_nx, &
      patch%coarse_ny, patch%coarse_nz))
    allocate(candidate_coarse_temperature( &
      patch%coarse_nx, patch%coarse_ny, patch%coarse_nz))
    do k = 1, patch%coarse_nz
      do j = 1, patch%coarse_ny
        do i = 1, patch%coarse_nx
          read(unit, *, iostat=status) &
            candidate_coarse_state(:, i, j, k), &
            candidate_coarse_temperature(i, j, k)
          if (status /= 0) go to 900
        end do
      end do
    end do
    if (.not. all(ieee_is_finite(candidate_coarse_state)) .or. &
        .not. all(ieee_is_finite(candidate_coarse_temperature))) go to 900
    if (minval(candidate_coarse_temperature) <= 0.0_dp) go to 900
    allocate(recovered_coarse_temperature( &
      patch%coarse_nx, patch%coarse_ny, patch%coarse_nz))
    call recover_reactive_temperatures_3d( &
      species, candidate_coarse_state, candidate_coarse_temperature, &
      patch%coarse_nx, patch%coarse_ny, patch%coarse_nz, &
      recovered_coarse_temperature, local_ok)
    if (.not. local_ok) go to 900
    temperature_consistent = all(checkpoint_real_matches( &
      candidate_coarse_temperature, recovered_coarse_temperature))
    if (.not. temperature_consistent) go to 900
    candidate_coarse_temperature = recovered_coarse_temperature

    read(unit, *, iostat=status) stored_fine_shape
    if (status /= 0) go to 900
    if (any(stored_fine_shape /= &
        [patch%fine_nx(), patch%fine_ny(), patch%fine_nz()])) go to 900
    allocate(candidate_fine_state(nvar, patch%fine_nx(), &
      patch%fine_ny(), patch%fine_nz()))
    allocate(candidate_fine_temperature( &
      patch%fine_nx(), patch%fine_ny(), patch%fine_nz()))
    do k = 1, patch%fine_nz()
      do j = 1, patch%fine_ny()
        do i = 1, patch%fine_nx()
          read(unit, *, iostat=status) &
            candidate_fine_state(:, i, j, k), &
            candidate_fine_temperature(i, j, k)
          if (status /= 0) go to 900
        end do
      end do
    end do
    if (.not. all(ieee_is_finite(candidate_fine_state)) .or. &
        .not. all(ieee_is_finite(candidate_fine_temperature))) go to 900
    if (minval(candidate_fine_temperature) <= 0.0_dp) go to 900
    allocate(recovered_fine_temperature( &
      patch%fine_nx(), patch%fine_ny(), patch%fine_nz()))
    call recover_reactive_temperatures_3d( &
      species, candidate_fine_state, candidate_fine_temperature, &
      patch%fine_nx(), patch%fine_ny(), patch%fine_nz(), &
      recovered_fine_temperature, local_ok)
    if (.not. local_ok) go to 900
    temperature_consistent = all(checkpoint_real_matches( &
      candidate_fine_temperature, recovered_fine_temperature))
    if (.not. temperature_consistent) go to 900
    candidate_fine_temperature = recovered_fine_temperature
    call checkpoint_levels_synchronized( &
      patch, candidate_coarse_state, candidate_fine_state, synchronized)
    if (.not. synchronized) go to 900
    read(unit, '(a)', iostat=status) end_marker
    if (status /= 0) go to 900
    if (trim(end_marker) /= "END_CHECKPOINT") go to 900
    read(unit, '(a)', iostat=status) magic
    if (.not. is_iostat_end(status)) go to 900
    close(unit, iostat=status)
    if (status /= 0) then
      message = "Failed to close 3D AMR checkpoint"
      return
    end if

    coarse_state = candidate_coarse_state
    coarse_temperature = candidate_coarse_temperature
    fine_state = candidate_fine_state
    fine_temperature = candidate_fine_temperature
    time = stored_time
    steps = stored_steps
    initial_integrals = candidate_integrals
    maximum_reflux = stored_maximum_reflux
    if (flow_config%transport_enabled) then
      maximum_transport_diffusivity = &
        candidate_maximum_transport_diffusivity
      minimum_transport_theta = candidate_minimum_transport_theta
    end if
    ok = .true.
    return

900 continue
    close(unit)
    if (len_trim(message) == 0) then
      message = "3D AMR checkpoint is truncated, invalid, or incompatible"
    end if
  end subroutine read_amr_reactive_3d_checkpoint

  subroutine write_transport_checkpoint_context( &
      unit, transport, flow_config, maximum_diffusivity, minimum_theta, status)
    integer, intent(in) :: unit
    type(gas_transport_species), intent(in) :: transport(:)
    type(reactive_3d_config), intent(in) :: flow_config
    real(dp), intent(in) :: maximum_diffusivity, minimum_theta
    integer, intent(out) :: status

    integer :: species_index

    status = 0
    write(unit, '(a)', iostat=status) "TRANSPORT_CONTEXT"
    if (status /= 0) return
    write(unit, '(a)', iostat=status) amr_reactive_3d_transport_operator
    if (status /= 0) return
    write(unit, '(a)', iostat=status) transport_parameter_convention
    if (status /= 0) return
    write(unit, '(a)', iostat=status) transport_checkpoint_phase
    if (status /= 0) return
    write(unit, '(i0)', iostat=status) size(transport)
    if (status /= 0) return
    do species_index = 1, size(transport)
      write(unit, '(a)', iostat=status) trim(transport(species_index)%name)
      if (status /= 0) return
      write(unit, '(i0)', iostat=status) transport(species_index)%geometry
      if (status /= 0) return
      write(unit, '(*(es27.18e3,1x))', iostat=status) &
        transport(species_index)%well_depth, &
        transport(species_index)%diameter, &
        transport(species_index)%dipole, &
        transport(species_index)%polarizability, &
        transport(species_index)%rotational_relaxation
      if (status /= 0) return
    end do
    write(unit, '(a)', iostat=status) "TRANSPORT_CONTROLS"
    if (status /= 0) return
    write(unit, '(*(i0,1x))', iostat=status) &
      merge(1, 0, flow_config%viscosity_enabled), &
      merge(1, 0, flow_config%thermal_conduction_enabled), &
      merge(1, 0, flow_config%species_diffusion_enabled), &
      merge(1, 0, flow_config%barodiffusion_enabled)
    if (status /= 0) return
    write(unit, '(es27.18e3)', iostat=status) flow_config%transport_cfl
    if (status /= 0) return
    write(unit, '(a)', iostat=status) "TRANSPORT_DIAGNOSTICS"
    if (status /= 0) return
    write(unit, '(*(es27.18e3,1x))', iostat=status) &
      maximum_diffusivity, minimum_theta
    if (status /= 0) return
    write(unit, '(a)', iostat=status) "END_TRANSPORT_CONTEXT"
  end subroutine write_transport_checkpoint_context

  subroutine read_transport_checkpoint_context( &
      unit, species, transport, flow_config, maximum_diffusivity, &
      minimum_theta, ok, message)
    integer, intent(in) :: unit
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(reactive_3d_config), intent(in) :: flow_config
    real(dp), intent(out) :: maximum_diffusivity, minimum_theta
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    real(dp) :: stored_values(5), expected_values(5)
    real(dp) :: stored_transport_cfl, stored_diagnostics(2)
    character(len=1024) :: marker, stored_name
    integer :: stored_count, stored_geometry, stored_controls(4)
    integer :: expected_controls(4), species_index, status

    ok = .false.
    message = ""
    maximum_diffusivity = 0.0_dp
    minimum_theta = 1.0_dp
    read(unit, '(a)', iostat=status) marker
    if (status /= 0) then
      message = "Transport 3D AMR checkpoint context is invalid"
      return
    end if
    if (trim(marker) /= "TRANSPORT_CONTEXT") then
      message = "Transport 3D AMR checkpoint context is invalid"
      return
    end if
    read(unit, '(a)', iostat=status) marker
    if (status /= 0) then
      message = "Transport 3D AMR checkpoint operator mismatch"
      return
    end if
    if (trim(marker) /= amr_reactive_3d_transport_operator) then
      message = "Transport 3D AMR checkpoint operator mismatch"
      return
    end if
    read(unit, '(a)', iostat=status) marker
    if (status /= 0) then
      message = "Transport 3D AMR checkpoint convention mismatch"
      return
    end if
    if (trim(marker) /= transport_parameter_convention) then
      message = "Transport 3D AMR checkpoint convention mismatch"
      return
    end if
    read(unit, '(a)', iostat=status) marker
    if (status /= 0) then
      message = "Transport 3D AMR checkpoint phase mismatch"
      return
    end if
    if (trim(marker) /= transport_checkpoint_phase) then
      message = "Transport 3D AMR checkpoint phase mismatch"
      return
    end if
    read(unit, *, iostat=status) stored_count
    if (status /= 0) then
      message = "Transport 3D AMR checkpoint database mismatch"
      return
    end if
    if (stored_count /= size(transport) .or. &
        stored_count /= size(species)) then
      message = "Transport 3D AMR checkpoint database mismatch"
      return
    end if
    do species_index = 1, size(transport)
      read(unit, '(a)', iostat=status) stored_name
      if (status /= 0) then
        message = "Transport 3D AMR checkpoint database mismatch"
        return
      end if
      if (trim(stored_name) /= trim(transport(species_index)%name) .or. &
          trim(stored_name) /= &
          trim(species(species_index)%name)) then
        message = "Transport 3D AMR checkpoint database mismatch"
        return
      end if
      read(unit, *, iostat=status) stored_geometry
      if (status /= 0) then
        message = "Transport 3D AMR checkpoint database mismatch"
        return
      end if
      if (stored_geometry /= transport(species_index)%geometry) then
        message = "Transport 3D AMR checkpoint database mismatch"
        return
      end if
      read(unit, *, iostat=status) stored_values
      if (status /= 0) then
        message = "Transport 3D AMR checkpoint database mismatch"
        return
      end if
      expected_values = [ &
        transport(species_index)%well_depth, &
        transport(species_index)%diameter, &
        transport(species_index)%dipole, &
        transport(species_index)%polarizability, &
        transport(species_index)%rotational_relaxation]
      if (.not. all(checkpoint_real_matches( &
          stored_values, expected_values))) then
        message = "Transport 3D AMR checkpoint database mismatch"
        return
      end if
    end do
    read(unit, '(a)', iostat=status) marker
    if (status /= 0) then
      message = "Transport 3D AMR checkpoint controls are invalid"
      return
    end if
    if (trim(marker) /= "TRANSPORT_CONTROLS") then
      message = "Transport 3D AMR checkpoint controls are invalid"
      return
    end if
    read(unit, *, iostat=status) stored_controls
    if (status /= 0) then
      message = "Transport 3D AMR checkpoint controls mismatch"
      return
    end if
    expected_controls = [ &
      merge(1, 0, flow_config%viscosity_enabled), &
      merge(1, 0, flow_config%thermal_conduction_enabled), &
      merge(1, 0, flow_config%species_diffusion_enabled), &
      merge(1, 0, flow_config%barodiffusion_enabled)]
    if (any(stored_controls /= expected_controls)) then
      message = "Transport 3D AMR checkpoint controls mismatch"
      return
    end if
    read(unit, *, iostat=status) stored_transport_cfl
    if (status /= 0) then
      message = "Transport 3D AMR checkpoint CFL mismatch"
      return
    end if
    if (.not. checkpoint_real_matches( &
        stored_transport_cfl, flow_config%transport_cfl)) then
      message = "Transport 3D AMR checkpoint CFL mismatch"
      return
    end if
    read(unit, '(a)', iostat=status) marker
    if (status /= 0) then
      message = "Transport 3D AMR checkpoint diagnostics are invalid"
      return
    end if
    if (trim(marker) /= "TRANSPORT_DIAGNOSTICS") then
      message = "Transport 3D AMR checkpoint diagnostics are invalid"
      return
    end if
    read(unit, *, iostat=status) stored_diagnostics
    if (status /= 0) then
      message = "Transport 3D AMR checkpoint diagnostics are invalid"
      return
    end if
    if (.not. all(ieee_is_finite(stored_diagnostics))) then
      message = "Transport 3D AMR checkpoint diagnostics are invalid"
      return
    end if
    if (stored_diagnostics(1) < 0.0_dp .or. &
        stored_diagnostics(2) < 0.0_dp .or. &
        stored_diagnostics(2) > 1.0_dp) then
      message = "Transport 3D AMR checkpoint diagnostics are invalid"
      return
    end if
    read(unit, '(a)', iostat=status) marker
    if (status /= 0) then
      message = "Transport 3D AMR checkpoint context is invalid"
      return
    end if
    if (trim(marker) /= "END_TRANSPORT_CONTEXT") then
      message = "Transport 3D AMR checkpoint context is invalid"
      return
    end if
    maximum_diffusivity = stored_diagnostics(1)
    minimum_theta = stored_diagnostics(2)
    ok = .true.
  end subroutine read_transport_checkpoint_context

  pure logical function valid_transport_checkpoint_context( &
      species, flow_config, transport, maximum_diffusivity, minimum_theta) &
      result(valid)
    type(nasa7_species), intent(in) :: species(:)
    type(reactive_3d_config), intent(in) :: flow_config
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: maximum_diffusivity, minimum_theta

    valid = flow_config%transport_enabled .and. &
      compatible_transport_database(species, transport)
    if (valid) valid = ieee_is_finite(flow_config%transport_cfl)
    if (valid) valid = flow_config%transport_cfl > 0.0_dp
    if (valid) valid = ieee_is_finite(maximum_diffusivity)
    if (valid) valid = maximum_diffusivity >= 0.0_dp
    if (valid) valid = ieee_is_finite(minimum_theta)
    if (valid) valid = minimum_theta >= 0.0_dp .and. &
      minimum_theta <= 1.0_dp
    if (valid) valid = .not. flow_config%barodiffusion_enabled .or. &
      flow_config%species_diffusion_enabled
  end function valid_transport_checkpoint_context

  subroutine write_selected_checkpoint_context( &
      unit, reactions, bundle_sha256, chemistry_integrator, &
      base_mole_fractions, relative_tolerance, absolute_tolerance, status)
    integer, intent(in) :: unit
    type(elementary_reaction), intent(in) :: reactions(:)
    character(len=*), intent(in) :: bundle_sha256, chemistry_integrator
    real(dp), intent(in) :: base_mole_fractions(:)
    real(dp), intent(in) :: relative_tolerance, absolute_tolerance
    integer, intent(out) :: status

    integer :: reaction_index, efficiency_size

    status = 0
    write(unit, '(a)', iostat=status) "SELECTED_CONTEXT"
    if (status /= 0) return
    write(unit, '(a)', iostat=status) trim(bundle_sha256)
    if (status /= 0) return
    write(unit, '(a)', iostat=status) trim(chemistry_integrator)
    if (status /= 0) return
    write(unit, '(i0)', iostat=status) size(base_mole_fractions)
    if (status /= 0) return
    write(unit, '(*(es27.18e3,1x))', iostat=status) base_mole_fractions
    if (status /= 0) return
    write(unit, '(a)', iostat=status) "REACTIONS"
    if (status /= 0) return
    write(unit, '(i0)', iostat=status) size(reactions)
    if (status /= 0) return
    do reaction_index = 1, size(reactions)
      efficiency_size = 0
      if (allocated(reactions(reaction_index)%third_body_efficiencies)) then
        efficiency_size = &
          size(reactions(reaction_index)%third_body_efficiencies)
      end if
      write(unit, '(a)', iostat=status) &
        trim(reactions(reaction_index)%equation)
      if (status /= 0) return
      write(unit, '(*(i0,1x))', iostat=status) &
        reactions(reaction_index)%kind, &
        merge(1, 0, reactions(reaction_index)%reversible), &
        merge(1, 0, &
          allocated(reactions(reaction_index)%third_body_efficiencies)), &
        efficiency_size, &
        size(reactions(reaction_index)%reactant_stoich), &
        size(reactions(reaction_index)%product_stoich)
      if (status /= 0) return
      write(unit, '(*(es27.18e3,1x))', iostat=status) &
        reactions(reaction_index)%reactant_stoich
      if (status /= 0) return
      write(unit, '(*(es27.18e3,1x))', iostat=status) &
        reactions(reaction_index)%product_stoich
      if (status /= 0) return
      write(unit, '(*(es27.18e3,1x))', iostat=status) &
        reactions(reaction_index)%forward_rate%pre_exponential, &
        reactions(reaction_index)%forward_rate%temperature_exponent, &
        reactions(reaction_index)%forward_rate%activation_energy, &
        reactions(reaction_index)%low_pressure_rate%pre_exponential, &
        reactions(reaction_index)%low_pressure_rate%temperature_exponent, &
        reactions(reaction_index)%low_pressure_rate%activation_energy, &
        reactions(reaction_index)%high_pressure_rate%pre_exponential, &
        reactions(reaction_index)%high_pressure_rate%temperature_exponent, &
        reactions(reaction_index)%high_pressure_rate%activation_energy
      if (status /= 0) return
      write(unit, '(i0)', iostat=status) &
        merge(1, 0, reactions(reaction_index)%troe%enabled)
      if (status /= 0) return
      write(unit, '(*(es27.18e3,1x))', iostat=status) &
        reactions(reaction_index)%troe%alpha, &
        reactions(reaction_index)%troe%temperature_3, &
        reactions(reaction_index)%troe%temperature_1, &
        reactions(reaction_index)%troe%temperature_2
      if (status /= 0) return
      if (allocated(reactions(reaction_index)%third_body_efficiencies)) then
        write(unit, '(*(es27.18e3,1x))', iostat=status) &
          reactions(reaction_index)%third_body_efficiencies
        if (status /= 0) return
      end if
    end do
    write(unit, '(a)', iostat=status) "CHEMISTRY_CONTROLS"
    if (status /= 0) return
    write(unit, '(*(es27.18e3,1x))', iostat=status) &
      relative_tolerance, absolute_tolerance
    if (status /= 0) return
    write(unit, '(a)', iostat=status) "END_SELECTED_CONTEXT"
  end subroutine write_selected_checkpoint_context

  subroutine read_selected_checkpoint_context( &
      unit, reactions, bundle_sha256, chemistry_integrator, &
      base_mole_fractions, relative_tolerance, absolute_tolerance, &
      ok, message)
    integer, intent(in) :: unit
    type(elementary_reaction), intent(in) :: reactions(:)
    character(len=*), intent(in) :: bundle_sha256, chemistry_integrator
    real(dp), intent(in) :: base_mole_fractions(:)
    real(dp), intent(in) :: relative_tolerance, absolute_tolerance
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    real(dp), allocatable :: stored_composition(:)
    real(dp), allocatable :: stored_reactants(:), stored_products(:)
    real(dp), allocatable :: stored_efficiencies(:)
    real(dp) :: stored_rates(9), expected_rates(9)
    real(dp) :: stored_troe(4), expected_troe(4)
    real(dp) :: stored_controls(2), expected_controls(2)
    character(len=1024) :: marker, stored_bundle, stored_integrator
    character(len=1024) :: stored_equation
    integer :: status, stored_composition_size, stored_reaction_count
    integer :: stored_header(6), expected_header(6)
    integer :: stored_troe_enabled, reaction_index, efficiency_size

    ok = .false.
    message = ""
    read(unit, '(a)', iostat=status) marker
    if (status /= 0) then
      message = "Selected 3D AMR checkpoint context is invalid"
      return
    end if
    if (trim(marker) /= "SELECTED_CONTEXT") then
      message = "Selected 3D AMR checkpoint context is invalid"
      return
    end if
    read(unit, '(a)', iostat=status) stored_bundle
    if (status /= 0) then
      message = "Selected 3D AMR checkpoint context is invalid"
      return
    end if
    if (.not. valid_bundle_sha256(stored_bundle)) then
      message = "Selected 3D AMR checkpoint context is invalid"
      return
    end if
    if (trim(stored_bundle) /= trim(bundle_sha256)) then
      message = "Selected 3D AMR checkpoint bundle SHA-256 mismatch"
      return
    end if
    read(unit, '(a)', iostat=status) stored_integrator
    if (status /= 0) then
      message = "Selected 3D AMR checkpoint context is invalid"
      return
    end if
    if (.not. valid_chemistry_integrator(stored_integrator)) then
      message = "Selected 3D AMR checkpoint context is invalid"
      return
    end if
    if (trim(stored_integrator) /= trim(chemistry_integrator)) then
      message = "Selected 3D AMR checkpoint chemistry integrator mismatch"
      return
    end if
    read(unit, *, iostat=status) stored_composition_size
    if (status /= 0) then
      message = "Selected 3D AMR checkpoint composition is invalid"
      return
    end if
    if (stored_composition_size /= size(base_mole_fractions)) then
      message = "Selected 3D AMR checkpoint composition is invalid"
      return
    end if
    allocate(stored_composition(stored_composition_size))
    read(unit, *, iostat=status) stored_composition
    if (status /= 0) then
      message = "Selected 3D AMR checkpoint composition is invalid"
      return
    end if
    if (.not. all(ieee_is_finite(stored_composition))) then
      message = "Selected 3D AMR checkpoint composition is invalid"
      return
    end if
    if (.not. all(checkpoint_real_matches( &
        stored_composition, base_mole_fractions))) then
      message = "Selected 3D AMR checkpoint composition mismatch"
      return
    end if
    read(unit, '(a)', iostat=status) marker
    if (status /= 0) then
      message = "Selected 3D AMR checkpoint mechanism is invalid"
      return
    end if
    if (trim(marker) /= "REACTIONS") then
      message = "Selected 3D AMR checkpoint mechanism is invalid"
      return
    end if
    read(unit, *, iostat=status) stored_reaction_count
    if (status /= 0) then
      message = "Selected 3D AMR checkpoint mechanism mismatch"
      return
    end if
    if (stored_reaction_count /= size(reactions)) then
      message = "Selected 3D AMR checkpoint mechanism mismatch"
      return
    end if
    do reaction_index = 1, size(reactions)
      efficiency_size = 0
      if (allocated(reactions(reaction_index)%third_body_efficiencies)) then
        efficiency_size = &
          size(reactions(reaction_index)%third_body_efficiencies)
      end if
      read(unit, '(a)', iostat=status) stored_equation
      if (status /= 0) then
        message = "Selected 3D AMR checkpoint mechanism mismatch"
        return
      end if
      if (trim(stored_equation) /= trim(reactions(reaction_index)%equation)) then
        message = "Selected 3D AMR checkpoint mechanism mismatch"
        return
      end if
      read(unit, *, iostat=status) stored_header
      if (status /= 0) then
        message = "Selected 3D AMR checkpoint mechanism mismatch"
        return
      end if
      expected_header = [ &
        reactions(reaction_index)%kind, &
        merge(1, 0, reactions(reaction_index)%reversible), &
        merge(1, 0, &
          allocated(reactions(reaction_index)%third_body_efficiencies)), &
        efficiency_size, &
        size(reactions(reaction_index)%reactant_stoich), &
        size(reactions(reaction_index)%product_stoich)]
      if (any(stored_header /= expected_header)) then
        message = "Selected 3D AMR checkpoint mechanism mismatch"
        return
      end if
      allocate(stored_reactants(expected_header(5)))
      allocate(stored_products(expected_header(6)))
      read(unit, *, iostat=status) stored_reactants
      if (status /= 0) then
        message = "Selected 3D AMR checkpoint mechanism is invalid"
        return
      end if
      if (.not. all(ieee_is_finite(stored_reactants))) then
        message = "Selected 3D AMR checkpoint mechanism is invalid"
        return
      end if
      read(unit, *, iostat=status) stored_products
      if (status /= 0) then
        message = "Selected 3D AMR checkpoint mechanism is invalid"
        return
      end if
      if (.not. all(ieee_is_finite(stored_products))) then
        message = "Selected 3D AMR checkpoint mechanism is invalid"
        return
      end if
      if (.not. all(checkpoint_real_matches( &
          stored_reactants, reactions(reaction_index)%reactant_stoich)) .or. &
          .not. all(checkpoint_real_matches( &
            stored_products, reactions(reaction_index)%product_stoich))) then
        message = "Selected 3D AMR checkpoint mechanism mismatch"
        return
      end if
      read(unit, *, iostat=status) stored_rates
      if (status /= 0) then
        message = "Selected 3D AMR checkpoint mechanism mismatch"
        return
      end if
      expected_rates = [ &
        reactions(reaction_index)%forward_rate%pre_exponential, &
        reactions(reaction_index)%forward_rate%temperature_exponent, &
        reactions(reaction_index)%forward_rate%activation_energy, &
        reactions(reaction_index)%low_pressure_rate%pre_exponential, &
        reactions(reaction_index)%low_pressure_rate%temperature_exponent, &
        reactions(reaction_index)%low_pressure_rate%activation_energy, &
        reactions(reaction_index)%high_pressure_rate%pre_exponential, &
        reactions(reaction_index)%high_pressure_rate%temperature_exponent, &
        reactions(reaction_index)%high_pressure_rate%activation_energy]
      if (.not. all(ieee_is_finite(stored_rates))) then
        message = "Selected 3D AMR checkpoint mechanism mismatch"
        return
      end if
      if (.not. all(checkpoint_real_matches(stored_rates, expected_rates))) then
        message = "Selected 3D AMR checkpoint mechanism mismatch"
        return
      end if
      read(unit, *, iostat=status) stored_troe_enabled
      if (status /= 0) then
        message = "Selected 3D AMR checkpoint mechanism mismatch"
        return
      end if
      if (stored_troe_enabled /= &
          merge(1, 0, reactions(reaction_index)%troe%enabled)) then
        message = "Selected 3D AMR checkpoint mechanism mismatch"
        return
      end if
      read(unit, *, iostat=status) stored_troe
      if (status /= 0) then
        message = "Selected 3D AMR checkpoint mechanism mismatch"
        return
      end if
      expected_troe = [ &
        reactions(reaction_index)%troe%alpha, &
        reactions(reaction_index)%troe%temperature_3, &
        reactions(reaction_index)%troe%temperature_1, &
        reactions(reaction_index)%troe%temperature_2]
      if (.not. all(ieee_is_finite(stored_troe))) then
        message = "Selected 3D AMR checkpoint mechanism mismatch"
        return
      end if
      if (.not. all(checkpoint_real_matches(stored_troe, expected_troe))) then
        message = "Selected 3D AMR checkpoint mechanism mismatch"
        return
      end if
      if (efficiency_size > 0) then
        allocate(stored_efficiencies(efficiency_size))
        read(unit, *, iostat=status) stored_efficiencies
        if (status /= 0) then
          message = "Selected 3D AMR checkpoint mechanism mismatch"
          return
        end if
        if (.not. all(ieee_is_finite(stored_efficiencies))) then
          message = "Selected 3D AMR checkpoint mechanism mismatch"
          return
        end if
        if (.not. all(checkpoint_real_matches( &
            stored_efficiencies, &
            reactions(reaction_index)%third_body_efficiencies))) then
          message = "Selected 3D AMR checkpoint mechanism mismatch"
          return
        end if
        deallocate(stored_efficiencies)
      end if
      deallocate(stored_reactants, stored_products)
    end do
    read(unit, '(a)', iostat=status) marker
    if (status /= 0) then
      message = "Selected 3D AMR checkpoint controls are invalid"
      return
    end if
    if (trim(marker) /= "CHEMISTRY_CONTROLS") then
      message = "Selected 3D AMR checkpoint controls are invalid"
      return
    end if
    read(unit, *, iostat=status) stored_controls
    if (status /= 0) then
      message = "Selected 3D AMR checkpoint chemistry controls mismatch"
      return
    end if
    expected_controls = [relative_tolerance, absolute_tolerance]
    if (.not. all(ieee_is_finite(stored_controls))) then
      message = "Selected 3D AMR checkpoint chemistry controls mismatch"
      return
    end if
    if (.not. all(checkpoint_real_matches( &
        stored_controls, expected_controls))) then
      message = "Selected 3D AMR checkpoint chemistry controls mismatch"
      return
    end if
    read(unit, '(a)', iostat=status) marker
    if (status /= 0) then
      message = "Selected 3D AMR checkpoint context is invalid"
      return
    end if
    if (trim(marker) /= "END_SELECTED_CONTEXT") then
      message = "Selected 3D AMR checkpoint context is invalid"
      return
    end if
    ok = .true.
  end subroutine read_selected_checkpoint_context

  logical function valid_selected_checkpoint_context( &
      nspecies, bundle_sha256, chemistry_integrator, &
      base_mole_fractions, reactions, relative_tolerance, &
      absolute_tolerance) result(valid)
    integer, intent(in) :: nspecies
    character(len=*), intent(in) :: bundle_sha256, chemistry_integrator
    real(dp), intent(in) :: base_mole_fractions(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: relative_tolerance, absolute_tolerance

    integer :: reaction_index

    valid = nspecies > 0
    if (valid) valid = valid_bundle_sha256(bundle_sha256)
    if (valid) valid = valid_chemistry_integrator(chemistry_integrator)
    if (valid) valid = size(base_mole_fractions) == nspecies
    if (valid) valid = all(ieee_is_finite(base_mole_fractions))
    if (valid) valid = minval(base_mole_fractions) >= 0.0_dp
    if (valid) valid = abs(sum(base_mole_fractions) - 1.0_dp) <= 5.0e-10_dp
    if (valid) valid = size(reactions) > 0
    if (valid) valid = all(ieee_is_finite( &
      [relative_tolerance, absolute_tolerance]))
    if (valid) valid = relative_tolerance > 0.0_dp .and. &
      absolute_tolerance > 0.0_dp
    if (.not. valid) return
    do reaction_index = 1, size(reactions)
      if (.not. valid_reaction_record( &
          reactions(reaction_index), nspecies)) then
        valid = .false.
        return
      end if
    end do
  end function valid_selected_checkpoint_context

  logical function valid_reaction_record( &
      reaction, nspecies) result(valid)
    type(elementary_reaction), intent(in) :: reaction
    integer, intent(in) :: nspecies

    valid = valid_elementary_reaction(reaction, nspecies)
    if (.not. valid) return
    valid = len_trim(reaction%equation) > 0
    if (.not. valid) return
    valid = valid_arrhenius_rate(reaction%forward_rate)
    if (.not. valid) return
    valid = valid_arrhenius_rate(reaction%low_pressure_rate)
    if (.not. valid) return
    valid = valid_arrhenius_rate(reaction%high_pressure_rate)
    if (.not. valid) return
    valid = valid_troe_parameters(reaction%troe)
    if (.not. valid) return
    valid = all(ieee_is_finite([ &
      reaction%troe%alpha, reaction%troe%temperature_3, &
      reaction%troe%temperature_1, reaction%troe%temperature_2]))
    if (.not. valid) return
    if (allocated(reaction%third_body_efficiencies)) then
      valid = size(reaction%third_body_efficiencies) == nspecies
      if (.not. valid) return
      valid = all(ieee_is_finite( &
        reaction%third_body_efficiencies))
      if (valid) valid = &
        all(reaction%third_body_efficiencies >= 0.0_dp)
    end if
  end function valid_reaction_record

  pure logical function valid_bundle_sha256(value) result(valid)
    character(len=*), intent(in) :: value

    integer :: index

    valid = len_trim(value) == bundle_sha256_length
    if (.not. valid) return
    do index = 1, bundle_sha256_length
      select case (value(index:index))
      case ('0':'9', 'a':'f')
      case default
        valid = .false.
        return
      end select
    end do
  end function valid_bundle_sha256

  pure logical function valid_chemistry_integrator(value) result(valid)
    character(len=*), intent(in) :: value

    valid = trim(value) == "explicit" .or. trim(value) == "implicit"
  end function valid_chemistry_integrator

  pure logical function checkpoint_configuration_matches( &
      flow_config, amr_config, patch, selected_context) result(matches)
    type(reactive_3d_config), intent(in) :: flow_config
    type(amr_reactive_3d_config), intent(in) :: amr_config
    type(amr_patch_3d), intent(in) :: patch
    logical, intent(in) :: selected_context

    logical :: physics_matches

    if (selected_context) then
      physics_matches = trim(flow_config%thermo_model) == "selected"
    else
      physics_matches = .not. flow_config%chemistry_enabled
    end if

    matches = patch%is_strictly_interior() .and. &
      patch%coarse_nx == flow_config%nx .and. &
      patch%coarse_ny == flow_config%ny .and. &
      patch%coarse_nz == flow_config%nz .and. &
      patch%coarse_i_lower == amr_config%coarse_i_lower .and. &
      patch%coarse_i_upper == amr_config%coarse_i_upper .and. &
      patch%coarse_j_lower == amr_config%coarse_j_lower .and. &
      patch%coarse_j_upper == amr_config%coarse_j_upper .and. &
      patch%coarse_k_lower == amr_config%coarse_k_lower .and. &
      patch%coarse_k_upper == amr_config%coarse_k_upper .and. &
      patch%refinement_ratio == amr_config%refinement_ratio .and. &
      trim(flow_config%boundary_condition) == "periodic" .and. &
      (trim(flow_config%reconstruction) == "pcm" .or. &
       trim(flow_config%reconstruction) == "characteristic_plm") .and. &
      (trim(flow_config%limiter) == "minmod" .or. &
       trim(flow_config%limiter) == "mc") .and. physics_matches
  end function checkpoint_configuration_matches

  pure logical function valid_checkpoint_shapes( &
      species, patch, coarse_state, coarse_temperature, &
      fine_state, fine_temperature) result(valid)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: coarse_state(:, :, :, :)
    real(dp), intent(in) :: coarse_temperature(:, :, :)
    real(dp), intent(in) :: fine_state(:, :, :, :)
    real(dp), intent(in) :: fine_temperature(:, :, :)

    integer :: nvar

    nvar = reactive_nvar(size(species))
    valid = patch%is_strictly_interior() .and. &
      size(coarse_state, 1) == nvar .and. &
      size(coarse_state, 2) == patch%coarse_nx .and. &
      size(coarse_state, 3) == patch%coarse_ny .and. &
      size(coarse_state, 4) == patch%coarse_nz .and. &
      all(shape(coarse_temperature) == &
        [patch%coarse_nx, patch%coarse_ny, patch%coarse_nz]) .and. &
      size(fine_state, 1) == nvar .and. &
      size(fine_state, 2) == patch%fine_nx() .and. &
      size(fine_state, 3) == patch%fine_ny() .and. &
      size(fine_state, 4) == patch%fine_nz() .and. &
      all(shape(fine_temperature) == &
        [patch%fine_nx(), patch%fine_ny(), patch%fine_nz()])
  end function valid_checkpoint_shapes

  subroutine checkpoint_levels_synchronized( &
      patch, coarse_state, fine_state, synchronized)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(in) :: coarse_state(:, :, :, :)
    real(dp), intent(in) :: fine_state(:, :, :, :)
    logical, intent(out) :: synchronized

    real(dp), allocatable :: restricted(:, :, :, :)
    real(dp) :: scale, error
    logical :: local_ok
    integer :: covered_nx, covered_ny, covered_nz

    synchronized = .false.
    covered_nx = patch%coarse_i_upper - patch%coarse_i_lower + 1
    covered_ny = patch%coarse_j_upper - patch%coarse_j_lower + 1
    covered_nz = patch%coarse_k_upper - patch%coarse_k_lower + 1
    allocate(restricted(size(coarse_state, 1), &
      covered_nx, covered_ny, covered_nz))
    call restrict_average_3d(fine_state, patch, restricted, local_ok)
    if (.not. local_ok) return
    scale = max(1.0_dp, maxval(abs(restricted)))
    error = maxval(abs(restricted - coarse_state(:, &
      patch%coarse_i_lower:patch%coarse_i_upper, &
      patch%coarse_j_lower:patch%coarse_j_upper, &
      patch%coarse_k_lower:patch%coarse_k_upper))) / scale
    synchronized = ieee_is_finite(error)
    if (synchronized) synchronized = &
      error <= 64.0_dp * epsilon(1.0_dp)
  end subroutine checkpoint_levels_synchronized

  subroutine checkpoint_temperature_is_consistent( &
      species, state, temperature, consistent)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    logical, intent(out) :: consistent

    real(dp), allocatable :: recovered_temperature(:, :, :)
    logical :: local_ok

    consistent = .false.
    allocate(recovered_temperature( &
      size(state, 2), size(state, 3), size(state, 4)))
    call recover_reactive_temperatures_3d( &
      species, state, temperature, size(state, 2), size(state, 3), &
      size(state, 4), recovered_temperature, local_ok)
    if (.not. local_ok) return
    consistent = all(checkpoint_real_matches( &
      temperature, recovered_temperature))
  end subroutine checkpoint_temperature_is_consistent

  pure elemental logical function checkpoint_real_matches( &
      stored, expected) result(matches)
    real(dp), intent(in) :: stored, expected

    real(dp) :: tolerance

    matches = ieee_is_finite(stored)
    if (.not. matches) return
    matches = ieee_is_finite(expected)
    if (.not. matches) return
    tolerance = 64.0_dp * epsilon(1.0_dp) * &
      max(1.0_dp, abs(stored), abs(expected))
    matches = abs(stored - expected) <= tolerance
  end function checkpoint_real_matches

end module amr_reactive_3d_checkpoint_mod
