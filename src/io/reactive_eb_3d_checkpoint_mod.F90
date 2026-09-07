module reactive_eb_3d_checkpoint_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use constants_mod, only: density_floor
  use state_indices_mod, only: irho, iet
  use nasa7_thermo_mod, only: nasa7_species, valid_nasa7_species
  use elementary_kinetics_mod, only: &
    elementary_reaction, valid_arrhenius_rate, valid_troe_parameters, &
    valid_elementary_reaction
  use gas_transport_mod, only: &
    gas_transport_species, valid_gas_transport_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_species_component, &
    reactive_conserved_to_primitive
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use eb_geometry_3d_mod, only: &
    eb_geometry_3d, eb_covered_cell_3d
  implicit none
  private

  character(len=*), parameter :: checkpoint_magic = &
    "PELEF_REACTIVE_EB_3D_CHECKPOINT"
  integer, parameter :: checkpoint_schema_fixed = 1
  integer, parameter :: checkpoint_schema_selected = 2
  integer, parameter :: bundle_sha256_length = 64
  integer, parameter :: n_config_reals = 19
  integer, parameter :: n_config_ints = 5
  integer, parameter :: n_config_flags = 7

  public :: write_reactive_eb_3d_checkpoint
  public :: read_reactive_eb_3d_checkpoint

contains

  subroutine write_reactive_eb_3d_checkpoint( &
      path, species, reactions, transport, config, geometry, state, &
      temperature, time, steps, initial_integrals, initial_l1_integrals, &
      initial_element_integrals, minimum_dt, maximum_transport_diffusivity, &
      minimum_transport_theta, ok, message, bundle_sha256, &
      chemistry_integrator, base_mole_fractions)
    character(len=*), intent(in) :: path
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(reactive_eb_3d_config), intent(in) :: config
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    real(dp), intent(in) :: time
    integer, intent(in) :: steps
    real(dp), intent(in) :: initial_integrals(:), initial_l1_integrals(:)
    real(dp), intent(in) :: initial_element_integrals(:)
    real(dp), intent(in) :: minimum_dt, maximum_transport_diffusivity
    real(dp), intent(in) :: minimum_transport_theta
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    character(len=*), intent(in), optional :: bundle_sha256
    character(len=*), intent(in), optional :: chemistry_integrator
    real(dp), intent(in), optional :: base_mole_fractions(:)

    integer :: unit, status, nvar, species_index, reaction_index
    integer :: efficiency_size
    integer :: transport_index, i, j, k, schema
    real(dp) :: diagnostic_values(3)
    logical :: selected_context

    ok = .false.
    message = ""
    unit = -1
    nvar = reactive_nvar(size(species))
    selected_context = present(bundle_sha256) .or. &
      present(chemistry_integrator) .or. present(base_mole_fractions)

    if (len_trim(path) == 0) then
      message = "Reactive EB 3D checkpoint path is empty"
      return
    end if
    if (.not. valid_species_array(species) .or. nvar <= 0) then
      message = "Invalid Reactive EB 3D checkpoint thermodynamics"
      return
    end if
    if (.not. valid_mechanism_array(reactions, size(species))) then
      message = "Invalid Reactive EB 3D checkpoint mechanism"
      return
    end if
    if (config%chemistry_enabled .and. size(reactions) == 0) then
      message = "Reactive EB 3D chemistry checkpoint has no reactions"
      return
    end if
    if (.not. valid_transport_array(transport, species)) then
      message = "Invalid Reactive EB 3D checkpoint transport"
      return
    end if
    if (config%transport_enabled .and. size(transport) /= size(species)) then
      message = "Reactive EB 3D transport checkpoint is incomplete"
      return
    end if
    if (.not. valid_checkpoint_config(config)) then
      message = "Invalid Reactive EB 3D checkpoint configuration"
      return
    end if
    if (.not. valid_checkpoint_geometry(geometry, config)) then
      message = "Invalid Reactive EB 3D checkpoint geometry"
      return
    end if
    if (.not. valid_checkpoint_shapes( &
          species, geometry, state, temperature, initial_integrals, &
          initial_l1_integrals, initial_element_integrals)) then
      message = "Invalid Reactive EB 3D checkpoint shapes"
      return
    end if
    if (.not. valid_checkpoint_state( &
          species, geometry, state, temperature)) then
      message = "Invalid Reactive EB 3D checkpoint state or temperature"
      return
    end if
    if (.not. valid_checkpoint_diagnostics( &
          config, time, steps, initial_integrals, initial_l1_integrals, &
          initial_element_integrals, minimum_dt, &
          maximum_transport_diffusivity, minimum_transport_theta)) then
      message = "Invalid Reactive EB 3D checkpoint diagnostics"
      return
    end if
    if (selected_context .and. .not. (present(bundle_sha256) .and. &
        present(chemistry_integrator) .and. &
        present(base_mole_fractions))) then
      message = "Reactive EB 3D selected checkpoint context is incomplete"
      return
    end if
    if (selected_context) then
      if (.not. valid_selected_checkpoint_context( &
          size(species), bundle_sha256, chemistry_integrator, &
          base_mole_fractions)) then
        message = "Reactive EB 3D selected checkpoint context is invalid"
        return
      end if
      schema = checkpoint_schema_selected
    else
      schema = checkpoint_schema_fixed
    end if

    open(newunit=unit, file=trim(path), status="replace", action="write", &
      form="formatted", iostat=status)
    if (status /= 0) then
      write(message, '(a,1x,a)') &
        "Could not create Reactive EB 3D checkpoint:", trim(path)
      return
    end if

    write(unit, '(a)', iostat=status) checkpoint_magic
    if (status /= 0) go to 900
    write(unit, '(*(i0,1x))', iostat=status) &
      schema, size(species), size(reactions), size(transport), nvar
    if (status /= 0) go to 900
    if (selected_context) then
      write(unit, '(a)', iostat=status) "SELECTED_CONTEXT"
      if (status /= 0) go to 900
      write(unit, '(a)', iostat=status) trim(bundle_sha256)
      if (status /= 0) go to 900
      write(unit, '(a)', iostat=status) trim(chemistry_integrator)
      if (status /= 0) go to 900
      write(unit, '(i0)', iostat=status) size(base_mole_fractions)
      if (status /= 0) go to 900
      write(unit, '(*(es27.18e3,1x))', iostat=status) base_mole_fractions
      if (status /= 0) go to 900
    end if

    write(unit, '(a)', iostat=status) "SPECIES"
    if (status /= 0) go to 900
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

    write(unit, '(a)', iostat=status) "REACTIONS"
    if (status /= 0) go to 900
    do reaction_index = 1, size(reactions)
      efficiency_size = 0
      if (allocated(reactions(reaction_index)%third_body_efficiencies)) &
        efficiency_size = &
          size(reactions(reaction_index)%third_body_efficiencies)
      write(unit, '(a)', iostat=status) trim(reactions(reaction_index)%equation)
      if (status /= 0) go to 900
      write(unit, '(*(i0,1x))', iostat=status) &
        reactions(reaction_index)%kind, &
        merge(1, 0, reactions(reaction_index)%reversible), &
        merge(1, 0, allocated(reactions(reaction_index)%third_body_efficiencies)), &
        efficiency_size, &
        size(reactions(reaction_index)%reactant_stoich), &
        size(reactions(reaction_index)%product_stoich)
      if (status /= 0) go to 900
      write(unit, '(*(es27.18e3,1x))', iostat=status) &
        reactions(reaction_index)%reactant_stoich
      if (status /= 0) go to 900
      write(unit, '(*(es27.18e3,1x))', iostat=status) &
        reactions(reaction_index)%product_stoich
      if (status /= 0) go to 900
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
      if (status /= 0) go to 900
      write(unit, '(*(i0,1x))', iostat=status) &
        merge(1, 0, reactions(reaction_index)%troe%enabled)
      if (status /= 0) go to 900
      write(unit, '(*(es27.18e3,1x))', iostat=status) &
        reactions(reaction_index)%troe%alpha, &
        reactions(reaction_index)%troe%temperature_3, &
        reactions(reaction_index)%troe%temperature_1, &
        reactions(reaction_index)%troe%temperature_2
      if (status /= 0) go to 900
      if (allocated(reactions(reaction_index)%third_body_efficiencies)) then
        write(unit, '(*(es27.18e3,1x))', iostat=status) &
          reactions(reaction_index)%third_body_efficiencies
        if (status /= 0) go to 900
      end if
    end do

    write(unit, '(a)', iostat=status) "TRANSPORT"
    if (status /= 0) go to 900
    do transport_index = 1, size(transport)
      write(unit, '(a)', iostat=status) trim(transport(transport_index)%name)
      if (status /= 0) go to 900
      write(unit, '(*(i0,1x))', iostat=status) &
        transport(transport_index)%geometry
      if (status /= 0) go to 900
      write(unit, '(*(es27.18e3,1x))', iostat=status) &
        transport(transport_index)%well_depth, &
        transport(transport_index)%diameter, &
        transport(transport_index)%dipole, &
        transport(transport_index)%polarizability, &
        transport(transport_index)%rotational_relaxation
      if (status /= 0) go to 900
    end do

    write(unit, '(a)', iostat=status) "CONFIG"
    if (status /= 0) go to 900
    write(unit, '(*(i0,1x))', iostat=status) &
      config%nx, config%ny, config%nz, config%maximum_steps, &
      config%checkpoint_interval_steps
    if (status /= 0) go to 900
    write(unit, '(*(es27.18e3,1x))', iostat=status) &
      config%x_lower, config%x_upper, config%y_lower, config%y_upper, &
      config%z_lower, config%z_upper, config%final_time, config%cfl, &
      config%plane_position, config%state_redist_target_volume_fraction, &
      config%chemistry_relative_tolerance, &
      config%chemistry_absolute_tolerance, config%transport_cfl, &
      config%regular_density, config%cut_density, config%initial_pressure, &
      config%initial_velocity_x, config%initial_velocity_y, &
      config%initial_velocity_z
    if (status /= 0) go to 900
    write(unit, '(a)', iostat=status) trim(config%plane_axis)
    if (status /= 0) go to 900
    write(unit, '(a)', iostat=status) trim(config%riemann_solver)
    if (status /= 0) go to 900
    write(unit, '(a)', iostat=status) trim(config%redistribution)
    if (status /= 0) go to 900
    write(unit, '(*(i0,1x))', iostat=status) &
      merge(1, 0, config%chemistry_enabled), &
      merge(1, 0, config%transport_enabled), &
      merge(1, 0, config%viscosity_enabled), &
      merge(1, 0, config%thermal_conduction_enabled), &
      merge(1, 0, config%species_diffusion_enabled), &
      merge(1, 0, config%barodiffusion_enabled), &
      merge(1, 0, config%stop_after_checkpoint)
    if (status /= 0) go to 900

    write(unit, '(a)', iostat=status) "GEOMETRY"
    if (status /= 0) go to 900
    write(unit, '(*(i0,1x))', iostat=status) &
      geometry%nx, geometry%ny, geometry%nz
    if (status /= 0) go to 900
    write(unit, '(*(es27.18e3,1x))', iostat=status) &
      geometry%x_lower, geometry%x_upper, geometry%y_lower, &
      geometry%y_upper, geometry%z_lower, geometry%z_upper, &
      geometry%dx, geometry%dy, geometry%dz
    if (status /= 0) go to 900
    call write_real_array_3d(unit, geometry%volume_fraction, status)
    if (status /= 0) go to 900
    call write_real_array_3d(unit, geometry%cell_centroid_x, status)
    if (status /= 0) go to 900
    call write_real_array_3d(unit, geometry%cell_centroid_y, status)
    if (status /= 0) go to 900
    call write_real_array_3d(unit, geometry%cell_centroid_z, status)
    if (status /= 0) go to 900
    call write_real_array_3d(unit, geometry%x_face_fraction, status)
    if (status /= 0) go to 900
    call write_real_array_3d(unit, geometry%y_face_fraction, status)
    if (status /= 0) go to 900
    call write_real_array_3d(unit, geometry%z_face_fraction, status)
    if (status /= 0) go to 900
    call write_real_array_3d(unit, geometry%x_face_centroid_y, status)
    if (status /= 0) go to 900
    call write_real_array_3d(unit, geometry%x_face_centroid_z, status)
    if (status /= 0) go to 900
    call write_real_array_3d(unit, geometry%y_face_centroid_x, status)
    if (status /= 0) go to 900
    call write_real_array_3d(unit, geometry%y_face_centroid_z, status)
    if (status /= 0) go to 900
    call write_real_array_3d(unit, geometry%z_face_centroid_x, status)
    if (status /= 0) go to 900
    call write_real_array_3d(unit, geometry%z_face_centroid_y, status)
    if (status /= 0) go to 900
    call write_real_array_3d(unit, geometry%boundary_area, status)
    if (status /= 0) go to 900
    call write_real_array_3d(unit, geometry%boundary_centroid_x, status)
    if (status /= 0) go to 900
    call write_real_array_3d(unit, geometry%boundary_centroid_y, status)
    if (status /= 0) go to 900
    call write_real_array_3d(unit, geometry%boundary_centroid_z, status)
    if (status /= 0) go to 900
    call write_real_array_3d(unit, geometry%boundary_normal_x, status)
    if (status /= 0) go to 900
    call write_real_array_3d(unit, geometry%boundary_normal_y, status)
    if (status /= 0) go to 900
    call write_real_array_3d(unit, geometry%boundary_normal_z, status)
    if (status /= 0) go to 900
    call write_real_array_3d(unit, geometry%boundary_normal_integral_x, status)
    if (status /= 0) go to 900
    call write_real_array_3d(unit, geometry%boundary_normal_integral_y, status)
    if (status /= 0) go to 900
    call write_real_array_3d(unit, geometry%boundary_normal_integral_z, status)
    if (status /= 0) go to 900
    call write_integer_array_3d(unit, geometry%cell_type, status)
    if (status /= 0) go to 900

    write(unit, '(a)', iostat=status) "DIAGNOSTICS"
    if (status /= 0) go to 900
    write(unit, '(es27.18e3,1x,i0)', iostat=status) time, steps
    if (status /= 0) go to 900
    write(unit, '(*(es27.18e3,1x))', iostat=status) initial_integrals
    if (status /= 0) go to 900
    write(unit, '(*(es27.18e3,1x))', iostat=status) initial_l1_integrals
    if (status /= 0) go to 900
    write(unit, '(*(es27.18e3,1x))', iostat=status) initial_element_integrals
    if (status /= 0) go to 900
    diagnostic_values = [ &
      minimum_dt, maximum_transport_diffusivity, minimum_transport_theta]
    write(unit, '(*(es27.18e3,1x))', iostat=status) diagnostic_values
    if (status /= 0) go to 900

    write(unit, '(a)', iostat=status) "STATE"
    if (status /= 0) go to 900
    write(unit, '(*(i0,1x))', iostat=status) &
      size(state, 2), size(state, 3), size(state, 4)
    if (status /= 0) go to 900
    do k = 1, size(state, 4)
      do j = 1, size(state, 3)
        do i = 1, size(state, 2)
          write(unit, '(*(es27.18e3,1x))', iostat=status) &
            state(:, i, j, k), temperature(i, j, k)
          if (status /= 0) go to 900
        end do
      end do
    end do

    write(unit, '(a)', iostat=status) "END_CHECKPOINT"
    if (status /= 0) go to 900
    close(unit, iostat=status)
    if (status /= 0) then
      message = "Failed to close Reactive EB 3D checkpoint"
      return
    end if
    ok = .true.
    return

900 continue
    close(unit, iostat=status)
    message = "Failed while writing Reactive EB 3D checkpoint"
  end subroutine write_reactive_eb_3d_checkpoint

  subroutine read_reactive_eb_3d_checkpoint( &
      path, species, reactions, transport, config, geometry, state, &
      temperature, time, steps, initial_integrals, initial_l1_integrals, &
      initial_element_integrals, minimum_dt, maximum_transport_diffusivity, &
      minimum_transport_theta, ok, message, bundle_sha256, &
      chemistry_integrator, base_mole_fractions)
    character(len=*), intent(in) :: path
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(reactive_eb_3d_config), intent(in) :: config
    type(eb_geometry_3d), intent(inout) :: geometry
    real(dp), intent(inout) :: state(:, :, :, :), temperature(:, :, :)
    real(dp), intent(inout) :: time
    integer, intent(inout) :: steps
    real(dp), intent(inout) :: initial_integrals(:), initial_l1_integrals(:)
    real(dp), intent(inout) :: initial_element_integrals(:)
    real(dp), intent(inout) :: minimum_dt, maximum_transport_diffusivity
    real(dp), intent(inout) :: minimum_transport_theta
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    character(len=*), intent(in), optional :: bundle_sha256
    character(len=*), intent(in), optional :: chemistry_integrator
    real(dp), intent(in), optional :: base_mole_fractions(:)

    type(eb_geometry_3d) :: candidate_geometry
    real(dp), allocatable :: candidate_state(:, :, :, :)
    real(dp), allocatable :: candidate_temperature(:, :, :)
    real(dp), allocatable :: candidate_integrals(:)
    real(dp), allocatable :: candidate_l1_integrals(:)
    real(dp) :: candidate_element_integrals(3)
    real(dp) :: stored_species_values(18), expected_species_values(18)
    real(dp) :: stored_rate_values(9), expected_rate_values(9)
    real(dp) :: stored_troe_values(4), expected_troe_values(4)
    real(dp) :: stored_transport_values(5), expected_transport_values(5)
    real(dp) :: stored_config_reals(n_config_reals)
    real(dp) :: stored_geometry_values(9)
    real(dp) :: stored_time, stored_minimum_dt
    real(dp) :: stored_maximum_transport_diffusivity
    real(dp) :: stored_minimum_transport_theta
    character(len=1024) :: magic, marker, stored_name, stored_equation
    character(len=1024) :: stored_axis, stored_solver, stored_redistribution
    character(len=1024) :: trailer
    integer :: stored_header(5), stored_reaction_header(6)
    integer :: stored_config_ints(n_config_ints)
    integer :: stored_config_flags(n_config_flags)
    integer :: stored_geometry_dims(3), stored_state_shape(3)
    integer :: stored_troe_enabled, stored_steps
    integer :: unit, status, nvar, species_index, reaction_index
    integer :: transport_index, i, j, k, l
    integer :: stored_kind, stored_reversible, stored_has_efficiency
    integer :: stored_efficiency_size, stored_reactant_size
    integer :: stored_product_size
    real(dp), allocatable :: stored_reactants(:), stored_products(:)
    real(dp), allocatable :: stored_efficiencies(:)
    real(dp), allocatable :: stored_base_mole_fractions(:)
    character(len=1024) :: stored_bundle_sha256
    character(len=1024) :: stored_chemistry_integrator
    integer :: expected_schema, stored_composition_size
    logical :: selected_context

    ok = .false.
    message = ""
    unit = -1
    nvar = reactive_nvar(size(species))
    selected_context = present(bundle_sha256) .or. &
      present(chemistry_integrator) .or. present(base_mole_fractions)

    if (len_trim(path) == 0) then
      message = "Reactive EB 3D restart path is empty"
      return
    end if
    if (.not. valid_species_array(species) .or. nvar <= 0) then
      message = "Invalid Reactive EB 3D restart target"
      return
    end if
    if (.not. valid_mechanism_array(reactions, size(species))) then
      message = "Invalid Reactive EB 3D restart target"
      return
    end if
    if (.not. valid_transport_array(transport, species)) then
      message = "Invalid Reactive EB 3D restart target"
      return
    end if
    if (config%chemistry_enabled .and. size(reactions) == 0) then
      message = "Invalid Reactive EB 3D restart target"
      return
    end if
    if (config%transport_enabled .and. &
        size(transport) /= size(species)) then
      message = "Invalid Reactive EB 3D restart target"
      return
    end if
    if (.not. valid_checkpoint_config(config)) then
      message = "Invalid Reactive EB 3D restart target"
      return
    end if
    if (.not. valid_checkpoint_geometry(geometry, config)) then
      message = "Invalid Reactive EB 3D restart target"
      return
    end if
    if (.not. valid_checkpoint_shapes( &
        species, geometry, state, temperature, initial_integrals, &
        initial_l1_integrals, initial_element_integrals)) then
      message = "Invalid Reactive EB 3D restart target"
      return
    end if
    if (selected_context .and. .not. (present(bundle_sha256) .and. &
        present(chemistry_integrator) .and. &
        present(base_mole_fractions))) then
      message = "Reactive EB 3D selected restart context is incomplete"
      return
    end if
    if (selected_context) then
      if (.not. valid_selected_checkpoint_context( &
          size(species), bundle_sha256, chemistry_integrator, &
          base_mole_fractions)) then
        message = "Reactive EB 3D selected restart context is invalid"
        return
      end if
      expected_schema = checkpoint_schema_selected
    else
      expected_schema = checkpoint_schema_fixed
    end if

    open(newunit=unit, file=trim(path), status="old", action="read", &
      form="formatted", iostat=status)
    if (status /= 0) then
      write(message, '(a,1x,a)') &
        "Could not open Reactive EB 3D checkpoint:", trim(path)
      return
    end if

    read(unit, '(a)', iostat=status) magic
    if (status /= 0) go to 900
    if (trim(magic) /= checkpoint_magic) go to 900
    read(unit, *, iostat=status) stored_header
    if (status /= 0) go to 900
    if (stored_header(1) /= expected_schema) then
      if (selected_context) then
        message = &
          "Reactive EB 3D checkpoint lacks selected mechanism context"
      else
        message = "Reactive EB 3D checkpoint is not a fixed-runtime schema"
      end if
      go to 900
    end if
    if (any(stored_header(2:5) /= &
        [size(species), size(reactions), size(transport), nvar])) go to 900

    if (selected_context) then
      read(unit, '(a)', iostat=status) marker
      if (status /= 0) go to 900
      if (trim(marker) /= "SELECTED_CONTEXT") then
        message = "Reactive EB 3D selected checkpoint context is invalid"
        go to 900
      end if
      read(unit, '(a)', iostat=status) stored_bundle_sha256
      if (status /= 0) go to 900
      if (.not. valid_bundle_sha256(stored_bundle_sha256)) then
        message = "Reactive EB 3D selected checkpoint context is invalid"
        go to 900
      end if
      if (trim(stored_bundle_sha256) /= trim(bundle_sha256)) then
        message = "Reactive EB 3D selected checkpoint bundle SHA-256 mismatch"
        go to 900
      end if
      read(unit, '(a)', iostat=status) stored_chemistry_integrator
      if (status /= 0) go to 900
      if (.not. valid_chemistry_integrator(stored_chemistry_integrator)) then
        message = "Reactive EB 3D selected checkpoint context is invalid"
        go to 900
      end if
      if (trim(stored_chemistry_integrator) /= &
          trim(chemistry_integrator)) then
        message = &
          "Reactive EB 3D selected checkpoint chemistry integrator mismatch"
        go to 900
      end if
      read(unit, *, iostat=status) stored_composition_size
      if (status /= 0) go to 900
      if (stored_composition_size /= size(species)) then
        message = "Reactive EB 3D selected checkpoint context is invalid"
        go to 900
      end if
      allocate(stored_base_mole_fractions(stored_composition_size))
      read(unit, *, iostat=status) stored_base_mole_fractions
      if (status /= 0) go to 900
      if (.not. valid_selected_composition( &
          stored_base_mole_fractions, size(species))) then
        message = "Reactive EB 3D selected checkpoint context is invalid"
        go to 900
      end if
      if (.not. all(checkpoint_real_matches( &
          stored_base_mole_fractions, base_mole_fractions))) then
        message = "Reactive EB 3D selected checkpoint composition mismatch"
        go to 900
      end if
    end if

    read(unit, '(a)', iostat=status) marker
    if (status /= 0) go to 900
    if (trim(marker) /= "SPECIES") go to 900
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
      if (.not. all(ieee_is_finite(stored_species_values)) .or. &
          .not. all(checkpoint_real_matches( &
            stored_species_values, expected_species_values))) go to 900
    end do

    read(unit, '(a)', iostat=status) marker
    if (status /= 0) go to 900
    if (trim(marker) /= "REACTIONS") go to 900
    do reaction_index = 1, size(reactions)
      read(unit, '(a)', iostat=status) stored_equation
      if (status /= 0) go to 900
      if (trim(stored_equation) /= trim(reactions(reaction_index)%equation)) go to 900
      read(unit, *, iostat=status) stored_reaction_header
      if (status /= 0) go to 900
      stored_kind = stored_reaction_header(1)
      stored_reversible = stored_reaction_header(2)
      stored_has_efficiency = stored_reaction_header(3)
      stored_efficiency_size = stored_reaction_header(4)
      stored_reactant_size = stored_reaction_header(5)
      stored_product_size = stored_reaction_header(6)
      if (stored_kind < 1 .or. stored_kind > 3 .or. &
          (stored_reversible /= 0 .and. stored_reversible /= 1) .or. &
          (stored_has_efficiency /= 0 .and. stored_has_efficiency /= 1) .or. &
          stored_reactant_size /= size(species) .or. &
          stored_product_size /= size(species) .or. &
          stored_efficiency_size < 0 .or. &
          stored_efficiency_size /= merge(size(species), 0, &
            stored_has_efficiency == 1)) go to 900
      if (stored_kind /= reactions(reaction_index)%kind .or. &
          stored_reversible /= merge(1, 0, reactions(reaction_index)%reversible) .or. &
          stored_has_efficiency /= merge(1, 0, &
            allocated(reactions(reaction_index)%third_body_efficiencies))) go to 900

      allocate(stored_reactants(size(species)), stored_products(size(species)))
      read(unit, *, iostat=status) stored_reactants
      if (status /= 0) go to 900
      if (.not. all(ieee_is_finite(stored_reactants))) go to 900
      read(unit, *, iostat=status) stored_products
      if (status /= 0) go to 900
      if (.not. all(ieee_is_finite(stored_products))) go to 900
      if (.not. all(checkpoint_real_matches(stored_reactants, &
          reactions(reaction_index)%reactant_stoich)) .or. &
          .not. all(checkpoint_real_matches(stored_products, &
            reactions(reaction_index)%product_stoich))) go to 900

      read(unit, *, iostat=status) stored_rate_values
      if (status /= 0) go to 900
      expected_rate_values = [ &
        reactions(reaction_index)%forward_rate%pre_exponential, &
        reactions(reaction_index)%forward_rate%temperature_exponent, &
        reactions(reaction_index)%forward_rate%activation_energy, &
        reactions(reaction_index)%low_pressure_rate%pre_exponential, &
        reactions(reaction_index)%low_pressure_rate%temperature_exponent, &
        reactions(reaction_index)%low_pressure_rate%activation_energy, &
        reactions(reaction_index)%high_pressure_rate%pre_exponential, &
        reactions(reaction_index)%high_pressure_rate%temperature_exponent, &
        reactions(reaction_index)%high_pressure_rate%activation_energy]
      if (.not. all(ieee_is_finite(stored_rate_values)) .or. &
          .not. all(checkpoint_real_matches( &
            stored_rate_values, expected_rate_values))) go to 900
      read(unit, *, iostat=status) stored_troe_enabled
      if (status /= 0) go to 900
      if ((stored_troe_enabled /= 0 .and. &
          stored_troe_enabled /= 1) .or. stored_troe_enabled /= &
          merge(1, 0, reactions(reaction_index)%troe%enabled)) go to 900
      read(unit, *, iostat=status) stored_troe_values
      if (status /= 0) go to 900
      expected_troe_values = [ &
        reactions(reaction_index)%troe%alpha, &
        reactions(reaction_index)%troe%temperature_3, &
        reactions(reaction_index)%troe%temperature_1, &
        reactions(reaction_index)%troe%temperature_2]
      if (.not. all(ieee_is_finite(stored_troe_values)) .or. &
          .not. all(checkpoint_real_matches( &
            stored_troe_values, expected_troe_values))) go to 900
      if (stored_has_efficiency == 1) then
        allocate(stored_efficiencies(size(species)))
        read(unit, *, iostat=status) stored_efficiencies
        if (status /= 0) go to 900
        if (.not. all(ieee_is_finite(stored_efficiencies)) .or. &
            .not. all(checkpoint_real_matches(stored_efficiencies, &
              reactions(reaction_index)%third_body_efficiencies))) go to 900
        deallocate(stored_efficiencies)
      end if
      deallocate(stored_reactants, stored_products)
    end do

    read(unit, '(a)', iostat=status) marker
    if (status /= 0) go to 900
    if (trim(marker) /= "TRANSPORT") go to 900
    do transport_index = 1, size(transport)
      read(unit, '(a)', iostat=status) stored_name
      if (status /= 0) go to 900
      if (trim(stored_name) /= trim(transport(transport_index)%name)) go to 900
      read(unit, *, iostat=status) stored_kind
      if (status /= 0) go to 900
      if (stored_kind /= transport(transport_index)%geometry) go to 900
      read(unit, *, iostat=status) stored_transport_values
      if (status /= 0) go to 900
      expected_transport_values = [ &
        transport(transport_index)%well_depth, &
        transport(transport_index)%diameter, &
        transport(transport_index)%dipole, &
        transport(transport_index)%polarizability, &
        transport(transport_index)%rotational_relaxation]
      if (.not. all(ieee_is_finite(stored_transport_values)) .or. &
          .not. all(checkpoint_real_matches( &
            stored_transport_values, expected_transport_values))) go to 900
    end do

    read(unit, '(a)', iostat=status) marker
    if (status /= 0) go to 900
    if (trim(marker) /= "CONFIG") go to 900
    read(unit, *, iostat=status) stored_config_ints
    if (status /= 0) go to 900
    read(unit, *, iostat=status) stored_config_reals
    if (status /= 0) go to 900
    if (.not. all(ieee_is_finite(stored_config_reals))) go to 900
    read(unit, '(a)', iostat=status) stored_axis
    if (status /= 0) go to 900
    read(unit, '(a)', iostat=status) stored_solver
    if (status /= 0) go to 900
    read(unit, '(a)', iostat=status) stored_redistribution
    if (status /= 0) go to 900
    read(unit, *, iostat=status) stored_config_flags
    if (status /= 0) go to 900
    if (.not. valid_flag_values(stored_config_flags)) go to 900
    if (.not. config_records_match( &
          config, stored_config_ints, stored_config_reals, stored_axis, &
          stored_solver, stored_redistribution, stored_config_flags)) go to 900

    read(unit, '(a)', iostat=status) marker
    if (status /= 0) go to 900
    if (trim(marker) /= "GEOMETRY") go to 900
    read(unit, *, iostat=status) stored_geometry_dims
    if (status /= 0) go to 900
    if (any(stored_geometry_dims /= &
        [geometry%nx, geometry%ny, geometry%nz])) go to 900
    candidate_geometry%nx = stored_geometry_dims(1)
    candidate_geometry%ny = stored_geometry_dims(2)
    candidate_geometry%nz = stored_geometry_dims(3)
    read(unit, *, iostat=status) stored_geometry_values
    if (status /= 0) go to 900
    if (.not. all(ieee_is_finite(stored_geometry_values))) go to 900
    candidate_geometry%x_lower = stored_geometry_values(1)
    candidate_geometry%x_upper = stored_geometry_values(2)
    candidate_geometry%y_lower = stored_geometry_values(3)
    candidate_geometry%y_upper = stored_geometry_values(4)
    candidate_geometry%z_lower = stored_geometry_values(5)
    candidate_geometry%z_upper = stored_geometry_values(6)
    candidate_geometry%dx = stored_geometry_values(7)
    candidate_geometry%dy = stored_geometry_values(8)
    candidate_geometry%dz = stored_geometry_values(9)
    call allocate_geometry_storage(candidate_geometry, status)
    if (status /= 0) go to 900
    call read_real_array_3d(unit, candidate_geometry%volume_fraction, status)
    if (status /= 0) go to 900
    call read_real_array_3d(unit, candidate_geometry%cell_centroid_x, status)
    if (status /= 0) go to 900
    call read_real_array_3d(unit, candidate_geometry%cell_centroid_y, status)
    if (status /= 0) go to 900
    call read_real_array_3d(unit, candidate_geometry%cell_centroid_z, status)
    if (status /= 0) go to 900
    call read_real_array_3d(unit, candidate_geometry%x_face_fraction, status)
    if (status /= 0) go to 900
    call read_real_array_3d(unit, candidate_geometry%y_face_fraction, status)
    if (status /= 0) go to 900
    call read_real_array_3d(unit, candidate_geometry%z_face_fraction, status)
    if (status /= 0) go to 900
    call read_real_array_3d(unit, candidate_geometry%x_face_centroid_y, status)
    if (status /= 0) go to 900
    call read_real_array_3d(unit, candidate_geometry%x_face_centroid_z, status)
    if (status /= 0) go to 900
    call read_real_array_3d(unit, candidate_geometry%y_face_centroid_x, status)
    if (status /= 0) go to 900
    call read_real_array_3d(unit, candidate_geometry%y_face_centroid_z, status)
    if (status /= 0) go to 900
    call read_real_array_3d(unit, candidate_geometry%z_face_centroid_x, status)
    if (status /= 0) go to 900
    call read_real_array_3d(unit, candidate_geometry%z_face_centroid_y, status)
    if (status /= 0) go to 900
    call read_real_array_3d(unit, candidate_geometry%boundary_area, status)
    if (status /= 0) go to 900
    call read_real_array_3d(unit, candidate_geometry%boundary_centroid_x, status)
    if (status /= 0) go to 900
    call read_real_array_3d(unit, candidate_geometry%boundary_centroid_y, status)
    if (status /= 0) go to 900
    call read_real_array_3d(unit, candidate_geometry%boundary_centroid_z, status)
    if (status /= 0) go to 900
    call read_real_array_3d(unit, candidate_geometry%boundary_normal_x, status)
    if (status /= 0) go to 900
    call read_real_array_3d(unit, candidate_geometry%boundary_normal_y, status)
    if (status /= 0) go to 900
    call read_real_array_3d(unit, candidate_geometry%boundary_normal_z, status)
    if (status /= 0) go to 900
    call read_real_array_3d(unit, candidate_geometry%boundary_normal_integral_x, status)
    if (status /= 0) go to 900
    call read_real_array_3d(unit, candidate_geometry%boundary_normal_integral_y, status)
    if (status /= 0) go to 900
    call read_real_array_3d(unit, candidate_geometry%boundary_normal_integral_z, status)
    if (status /= 0) go to 900
    call read_integer_array_3d(unit, candidate_geometry%cell_type, status)
    if (status /= 0) go to 900
    if (.not. candidate_geometry%is_valid()) go to 900
    if (.not. geometries_match(candidate_geometry, geometry)) go to 900

    read(unit, '(a)', iostat=status) marker
    if (status /= 0) go to 900
    if (trim(marker) /= "DIAGNOSTICS") go to 900
    read(unit, *, iostat=status) stored_time, stored_steps
    if (status /= 0) go to 900
    allocate(candidate_integrals(nvar), candidate_l1_integrals(nvar))
    read(unit, *, iostat=status) candidate_integrals
    if (status /= 0) go to 900
    read(unit, *, iostat=status) candidate_l1_integrals
    if (status /= 0) go to 900
    read(unit, *, iostat=status) candidate_element_integrals
    if (status /= 0) go to 900
    read(unit, *, iostat=status) stored_minimum_dt, &
      stored_maximum_transport_diffusivity, stored_minimum_transport_theta
    if (status /= 0) go to 900
    if (.not. all(ieee_is_finite(candidate_integrals)) .or. &
        .not. all(ieee_is_finite(candidate_l1_integrals)) .or. &
        .not. all(ieee_is_finite(candidate_element_integrals)) .or. &
        .not. all(ieee_is_finite([stored_time, stored_minimum_dt, &
          stored_maximum_transport_diffusivity, stored_minimum_transport_theta])) .or. &
        .not. valid_checkpoint_diagnostics( &
          config, stored_time, stored_steps, candidate_integrals, &
          candidate_l1_integrals, candidate_element_integrals, &
          stored_minimum_dt, stored_maximum_transport_diffusivity, &
          stored_minimum_transport_theta)) go to 900

    read(unit, '(a)', iostat=status) marker
    if (status /= 0) go to 900
    if (trim(marker) /= "STATE") go to 900
    read(unit, *, iostat=status) stored_state_shape
    if (status /= 0) go to 900
    if (any(stored_state_shape /= &
        [geometry%nx, geometry%ny, geometry%nz])) go to 900
    allocate(candidate_state(nvar, geometry%nx, geometry%ny, geometry%nz))
    allocate(candidate_temperature(geometry%nx, geometry%ny, geometry%nz))
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          read(unit, *, iostat=status) candidate_state(:, i, j, k), &
            candidate_temperature(i, j, k)
          if (status /= 0) go to 900
        end do
      end do
    end do
    if (.not. valid_checkpoint_state( &
          species, candidate_geometry, candidate_state, candidate_temperature)) go to 900

    read(unit, '(a)', iostat=status) marker
    if (status /= 0) go to 900
    if (trim(marker) /= "END_CHECKPOINT") go to 900
    read(unit, '(a)', iostat=status) trailer
    if (.not. is_iostat_end(status)) go to 900
    close(unit, iostat=status)
    if (status /= 0) then
      message = "Failed to close Reactive EB 3D checkpoint"
      return
    end if

    ! Nothing above this point aliases a public restart target.  Publish the
    ! complete candidate only after every record and EOS check has passed.
    geometry = candidate_geometry
    state = candidate_state
    temperature = candidate_temperature
    time = stored_time
    steps = stored_steps
    initial_integrals = candidate_integrals
    initial_l1_integrals = candidate_l1_integrals
    initial_element_integrals = candidate_element_integrals
    minimum_dt = stored_minimum_dt
    maximum_transport_diffusivity = stored_maximum_transport_diffusivity
    minimum_transport_theta = stored_minimum_transport_theta
    ok = .true.
    return

900 continue
    close(unit, iostat=l)
    if (len_trim(message) == 0) then
      message = &
        "Reactive EB 3D checkpoint is truncated, invalid, or incompatible"
    end if
  end subroutine read_reactive_eb_3d_checkpoint

  logical function valid_selected_checkpoint_context( &
      nspecies, bundle_sha256, chemistry_integrator, base_mole_fractions) &
      result(valid)
    integer, intent(in) :: nspecies
    character(len=*), intent(in) :: bundle_sha256, chemistry_integrator
    real(dp), intent(in) :: base_mole_fractions(:)

    valid = valid_bundle_sha256(bundle_sha256) .and. &
      valid_chemistry_integrator(chemistry_integrator) .and. &
      valid_selected_composition(base_mole_fractions, nspecies)
  end function valid_selected_checkpoint_context

  logical function valid_bundle_sha256(value) result(valid)
    character(len=*), intent(in) :: value
    integer :: i

    valid = len_trim(value) == bundle_sha256_length
    if (.not. valid) return
    do i = 1, bundle_sha256_length
      if (index("0123456789abcdef", value(i:i)) == 0) then
        valid = .false.
        return
      end if
    end do
  end function valid_bundle_sha256

  logical function valid_chemistry_integrator(value) result(valid)
    character(len=*), intent(in) :: value

    valid = trim(value) == "explicit" .or. trim(value) == "implicit"
  end function valid_chemistry_integrator

  logical function valid_selected_composition(values, nspecies) result(valid)
    real(dp), intent(in) :: values(:)
    integer, intent(in) :: nspecies

    valid = size(values) == nspecies .and. nspecies > 0
    if (.not. valid) return
    valid = all(ieee_is_finite(values)) .and. all(values >= 0.0_dp)
    if (.not. valid) return
    valid = checkpoint_real_matches(sum(values), 1.0_dp)
  end function valid_selected_composition

  logical function valid_species_array(species) result(valid)
    type(nasa7_species), intent(in) :: species(:)
    integer :: i, j

    valid = size(species) > 0
    if (.not. valid) return
    do i = 1, size(species)
      if (.not. valid_nasa7_species(species(i)) .or. &
          len_trim(species(i)%name) == 0) then
        valid = .false.
        return
      end if
      do j = 1, i - 1
        if (trim(species(i)%name) == trim(species(j)%name)) then
          valid = .false.
          return
        end if
      end do
    end do
  end function valid_species_array

  logical function valid_mechanism_array(reactions, nspecies) result(valid)
    type(elementary_reaction), intent(in) :: reactions(:)
    integer, intent(in) :: nspecies
    integer :: i

    valid = nspecies > 0
    if (.not. valid) return
    do i = 1, size(reactions)
      if (.not. valid_reaction_record(reactions(i), nspecies)) then
        valid = .false.
        return
      end if
    end do
  end function valid_mechanism_array

  logical function valid_reaction_record(reaction, nspecies) result(valid)
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
      valid = size(reaction%third_body_efficiencies) == nspecies .and. &
        all(ieee_is_finite(reaction%third_body_efficiencies)) .and. &
        all(reaction%third_body_efficiencies >= 0.0_dp)
    end if
  end function valid_reaction_record

  logical function valid_transport_array(transport, species) result(valid)
    type(gas_transport_species), intent(in) :: transport(:)
    type(nasa7_species), intent(in) :: species(:)
    integer :: i

    valid = size(transport) == 0 .or. size(transport) == size(species)
    if (.not. valid) return
    do i = 1, size(transport)
      valid = valid_gas_transport_species(transport(i)) .and. &
        trim(transport(i)%name) == trim(species(i)%name)
      if (.not. valid) return
    end do
  end function valid_transport_array

  logical function valid_checkpoint_config(config) result(valid)
    type(reactive_eb_3d_config), intent(in) :: config
    real(dp) :: axis_lower, axis_upper, spacing, coordinate, tolerance
    real(dp) :: cut_fraction

    valid = min(config%nx, config%ny, config%nz) >= 2 .and. &
      config%maximum_steps > 0 .and. &
      all(ieee_is_finite([config%x_lower, config%x_upper, &
        config%y_lower, config%y_upper, config%z_lower, config%z_upper, &
        config%final_time, config%cfl, config%plane_position, &
        config%state_redist_target_volume_fraction, &
        config%chemistry_relative_tolerance, &
        config%chemistry_absolute_tolerance, config%transport_cfl, &
        config%regular_density, config%cut_density, config%initial_pressure, &
        config%initial_velocity_x, config%initial_velocity_y, &
        config%initial_velocity_z])) .and. &
      config%x_upper > config%x_lower .and. &
      config%y_upper > config%y_lower .and. &
      config%z_upper > config%z_lower .and. config%final_time > 0.0_dp .and. &
      config%cfl > 0.0_dp .and. config%cfl <= 1.0_dp .and. &
      config%state_redist_target_volume_fraction > 0.0_dp .and. &
      config%state_redist_target_volume_fraction <= 1.0_dp .and. &
      config%chemistry_relative_tolerance > 0.0_dp .and. &
      config%chemistry_absolute_tolerance > 0.0_dp .and. &
      config%transport_cfl > 0.0_dp .and. config%transport_cfl <= 0.5_dp .and. &
      config%regular_density > 0.0_dp .and. config%cut_density > 0.0_dp .and. &
      config%initial_pressure > 0.0_dp .and. config%checkpoint_interval_steps >= 0
    if (.not. valid) return
    if (config%transport_enabled .and. &
        trim(config%redistribution) /= "state_redist") then
      valid = .false.
      return
    end if
    if (config%transport_enabled .and. &
        .not. (config%viscosity_enabled .or. &
          config%thermal_conduction_enabled .or. &
          config%species_diffusion_enabled)) then
      valid = .false.
      return
    end if
    if (config%transport_enabled .and. config%barodiffusion_enabled .and. &
        .not. config%species_diffusion_enabled) then
      valid = .false.
      return
    end if
    if (config%stop_after_checkpoint .and. config%checkpoint_interval_steps <= 0) then
      valid = .false.
      return
    end if
    if (trim(config%riemann_solver) /= "rusanov" .and. &
        trim(config%riemann_solver) /= "hllc" .and. &
        trim(config%riemann_solver) /= "pelec") then
      valid = .false.
      return
    end if
    if (trim(config%redistribution) /= "flux_redist" .and. &
        trim(config%redistribution) /= "state_redist") then
      valid = .false.
      return
    end if
    select case (trim(config%plane_axis))
    case ("x")
      axis_lower = config%x_lower
      axis_upper = config%x_upper
      spacing = (axis_upper - axis_lower) / real(config%nx, dp)
    case ("y")
      axis_lower = config%y_lower
      axis_upper = config%y_upper
      spacing = (axis_upper - axis_lower) / real(config%ny, dp)
    case ("z")
      axis_lower = config%z_lower
      axis_upper = config%z_upper
      spacing = (axis_upper - axis_lower) / real(config%nz, dp)
    case default
      valid = .false.
      return
    end select
    if (config%plane_position <= axis_lower .or. &
        config%plane_position >= axis_upper - spacing) then
      valid = .false.
      return
    end if
    coordinate = (config%plane_position - axis_lower) / spacing
    tolerance = 512.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(coordinate))
    if (abs(coordinate - anint(coordinate)) <= tolerance) then
      valid = .false.
      return
    end if
    cut_fraction = real(ceiling(coordinate), dp) - coordinate
    if (trim(config%redistribution) == "state_redist" .and. &
        config%state_redist_target_volume_fraction <= cut_fraction + tolerance) then
      valid = .false.
    end if
  end function valid_checkpoint_config

  logical function valid_checkpoint_geometry(geometry, config) result(valid)
    type(eb_geometry_3d), intent(in) :: geometry
    type(reactive_eb_3d_config), intent(in) :: config

    valid = geometry%is_valid() .and. geometry%nx == config%nx .and. &
      geometry%ny == config%ny .and. geometry%nz == config%nz .and. &
      all(checkpoint_real_matches([geometry%x_lower, geometry%x_upper, &
        geometry%y_lower, geometry%y_upper, geometry%z_lower, geometry%z_upper], &
        [config%x_lower, config%x_upper, config%y_lower, config%y_upper, &
          config%z_lower, config%z_upper]))
  end function valid_checkpoint_geometry

  logical function valid_checkpoint_shapes( &
      species, geometry, state, temperature, initial_integrals, &
      initial_l1_integrals, initial_element_integrals) result(valid)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    real(dp), intent(in) :: initial_integrals(:), initial_l1_integrals(:)
    real(dp), intent(in) :: initial_element_integrals(:)
    integer :: nvar

    nvar = reactive_nvar(size(species))
    valid = nvar > 0 .and. size(state, 1) == nvar .and. &
      size(state, 2) == geometry%nx .and. size(state, 3) == geometry%ny .and. &
      size(state, 4) == geometry%nz .and. all(shape(temperature) == &
        [geometry%nx, geometry%ny, geometry%nz]) .and. &
      size(initial_integrals) == nvar .and. &
      size(initial_l1_integrals) == nvar .and. &
      size(initial_element_integrals) == 3
  end function valid_checkpoint_shapes

  logical function valid_checkpoint_state( &
      species, geometry, state, temperature) result(valid)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    real(dp), allocatable :: primitive(:)
    real(dp) :: recovered_temperature, rho, species_total, closure
    real(dp) :: kinetic_density
    logical :: local_ok
    integer :: i, j, k, nvar

    valid = .false.
    nvar = reactive_nvar(size(species))
    if (size(state, 1) /= nvar .or. size(state, 2) /= geometry%nx .or. &
        size(state, 3) /= geometry%ny .or. size(state, 4) /= geometry%nz .or. &
        any(shape(temperature) /= [geometry%nx, geometry%ny, geometry%nz])) return
    if (.not. all(ieee_is_finite(state)) .or. &
        .not. all(ieee_is_finite(temperature)) .or. &
        minval(temperature) <= 0.0_dp) return
    allocate(primitive(5 + size(species)))
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          rho = state(irho, i, j, k)
          if (rho <= density_floor) return
          species_total = sum(state(reactive_species_component(1): &
            reactive_species_component(size(species)), i, j, k))
          closure = abs(species_total - rho) / max(rho, density_floor)
          if (species_total <= 0.0_dp .or. closure > 5.0e-10_dp .or. &
              any(state(reactive_species_component(1): &
                reactive_species_component(size(species)), i, j, k) < &
                -5.0e-10_dp * max(1.0_dp, rho))) return
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) then
            kinetic_density = 0.5_dp * &
              (state(2, i, j, k)**2 + state(3, i, j, k)**2 + &
               state(4, i, j, k)**2) / rho
            if (.not. ieee_is_finite(kinetic_density) .or. &
                .not. ieee_is_finite(state(iet, i, j, k))) return
            cycle
          end if
          call reactive_conserved_to_primitive( &
            species, state(:, i, j, k), temperature(i, j, k), primitive, &
            recovered_temperature, kinetic_density, local_ok)
          if (.not. local_ok .or. .not. checkpoint_real_matches( &
              temperature(i, j, k), recovered_temperature)) return
        end do
      end do
    end do
    valid = .true.
  end function valid_checkpoint_state

  logical function valid_checkpoint_diagnostics( &
      config, time, steps, initial_integrals, initial_l1_integrals, &
      initial_element_integrals, minimum_dt, maximum_transport_diffusivity, &
      minimum_transport_theta) result(valid)
    type(reactive_eb_3d_config), intent(in) :: config
    real(dp), intent(in) :: time, initial_integrals(:), initial_l1_integrals(:)
    real(dp), intent(in) :: initial_element_integrals(:)
    real(dp), intent(in) :: minimum_dt, maximum_transport_diffusivity
    real(dp), intent(in) :: minimum_transport_theta
    integer, intent(in) :: steps
    real(dp) :: time_tolerance

    time_tolerance = 64.0_dp * epsilon(1.0_dp) * &
      max(1.0_dp, abs(config%final_time))
    valid = size(initial_integrals) == size(initial_l1_integrals) .and. &
      size(initial_element_integrals) == 3 .and. &
      all(ieee_is_finite(initial_integrals)) .and. &
      all(ieee_is_finite(initial_l1_integrals)) .and. &
      all(ieee_is_finite(initial_element_integrals)) .and. &
      all(initial_l1_integrals >= 0.0_dp) .and. &
      all(initial_element_integrals >= 0.0_dp) .and. &
      ieee_is_finite(time) .and. time >= 0.0_dp .and. &
      time <= config%final_time + time_tolerance .and. steps >= 0 .and. &
      steps <= config%maximum_steps .and. &
      (steps > 0 .or. time <= tiny(1.0_dp)) .and. &
      (steps == 0 .or. time > 0.0_dp) .and. &
      ieee_is_finite(minimum_dt) .and. minimum_dt >= 0.0_dp .and. &
      (steps == 0 .or. minimum_dt > 0.0_dp) .and. &
      ieee_is_finite(maximum_transport_diffusivity) .and. &
      maximum_transport_diffusivity >= 0.0_dp .and. &
      ieee_is_finite(minimum_transport_theta) .and. &
      minimum_transport_theta >= 0.0_dp .and. minimum_transport_theta <= 1.0_dp
  ! The first expression above is intentionally replaced by the exact shape
  ! check below; it keeps this routine independent of the mechanism size.
    valid = valid .and. size(initial_integrals) > 0
  end function valid_checkpoint_diagnostics

  logical function config_records_match( &
      config, stored_ints, stored_reals, stored_axis, stored_solver, &
      stored_redistribution, stored_flags) result(matches)
    type(reactive_eb_3d_config), intent(in) :: config
    integer, intent(in) :: stored_ints(:), stored_flags(:)
    real(dp), intent(in) :: stored_reals(:)
    character(len=*), intent(in) :: stored_axis, stored_solver
    character(len=*), intent(in) :: stored_redistribution
    integer :: expected_ints(n_config_ints), expected_flags(n_config_flags)
    real(dp) :: expected_reals(n_config_reals)

    expected_ints = [config%nx, config%ny, config%nz, &
      config%maximum_steps, config%checkpoint_interval_steps]
    expected_reals = [ &
      config%x_lower, config%x_upper, config%y_lower, config%y_upper, &
      config%z_lower, config%z_upper, config%final_time, config%cfl, &
      config%plane_position, config%state_redist_target_volume_fraction, &
      config%chemistry_relative_tolerance, &
      config%chemistry_absolute_tolerance, config%transport_cfl, &
      config%regular_density, config%cut_density, config%initial_pressure, &
      config%initial_velocity_x, config%initial_velocity_y, &
      config%initial_velocity_z]
    expected_flags = [ &
      merge(1, 0, config%chemistry_enabled), &
      merge(1, 0, config%transport_enabled), &
      merge(1, 0, config%viscosity_enabled), &
      merge(1, 0, config%thermal_conduction_enabled), &
      merge(1, 0, config%species_diffusion_enabled), &
      merge(1, 0, config%barodiffusion_enabled), &
      merge(1, 0, config%stop_after_checkpoint)]
    matches = size(stored_ints) == n_config_ints .and. &
      size(stored_reals) == n_config_reals .and. &
      size(stored_flags) == n_config_flags .and. &
      stored_ints(1) == expected_ints(1) .and. &
      stored_ints(2) == expected_ints(2) .and. &
      stored_ints(3) == expected_ints(3) .and. &
      all(stored_flags(1:6) == expected_flags(1:6)) .and. &
      trim(stored_axis) == trim(config%plane_axis) .and. &
      trim(stored_solver) == trim(config%riemann_solver) .and. &
      trim(stored_redistribution) == trim(config%redistribution) .and. &
      all(checkpoint_real_matches(stored_reals(1:6), expected_reals(1:6))) .and. &
      all(checkpoint_real_matches(stored_reals(8:19), expected_reals(8:19)))
    ! The stored final_time and maximum_steps are provenance fields.  A
    ! restart may extend either bound; the stored time/step records are
    ! checked against the current bounds in valid_checkpoint_diagnostics.
    if (matches) then
      matches = stored_reals(7) > 0.0_dp .and. &
        ieee_is_finite(stored_reals(7)) .and. &
        stored_ints(4) > 0 .and. stored_ints(5) >= 0
    end if
  end function config_records_match

  logical function valid_flag_values(flags) result(valid)
    integer, intent(in) :: flags(:)
    valid = size(flags) == n_config_flags .and. &
      all((flags == 0) .or. (flags == 1))
  end function valid_flag_values

  logical function geometries_match(left, right) result(matches)
    type(eb_geometry_3d), intent(in) :: left, right

    matches = left%nx == right%nx .and. left%ny == right%ny .and. &
      left%nz == right%nz .and. &
      all(checkpoint_real_matches([left%x_lower, left%x_upper, left%y_lower, &
        left%y_upper, left%z_lower, left%z_upper, left%dx, left%dy, left%dz], &
        [right%x_lower, right%x_upper, right%y_lower, right%y_upper, &
          right%z_lower, right%z_upper, right%dx, right%dy, right%dz]))
    if (.not. matches) return
    matches = all(checkpoint_real_matches(left%volume_fraction, right%volume_fraction)) .and. &
      all(checkpoint_real_matches(left%cell_centroid_x, right%cell_centroid_x)) .and. &
      all(checkpoint_real_matches(left%cell_centroid_y, right%cell_centroid_y)) .and. &
      all(checkpoint_real_matches(left%cell_centroid_z, right%cell_centroid_z)) .and. &
      all(checkpoint_real_matches(left%x_face_fraction, right%x_face_fraction)) .and. &
      all(checkpoint_real_matches(left%y_face_fraction, right%y_face_fraction)) .and. &
      all(checkpoint_real_matches(left%z_face_fraction, right%z_face_fraction)) .and. &
      all(checkpoint_real_matches(left%x_face_centroid_y, right%x_face_centroid_y)) .and. &
      all(checkpoint_real_matches(left%x_face_centroid_z, right%x_face_centroid_z)) .and. &
      all(checkpoint_real_matches(left%y_face_centroid_x, right%y_face_centroid_x)) .and. &
      all(checkpoint_real_matches(left%y_face_centroid_z, right%y_face_centroid_z)) .and. &
      all(checkpoint_real_matches(left%z_face_centroid_x, right%z_face_centroid_x)) .and. &
      all(checkpoint_real_matches(left%z_face_centroid_y, right%z_face_centroid_y)) .and. &
      all(checkpoint_real_matches(left%boundary_area, right%boundary_area)) .and. &
      all(checkpoint_real_matches(left%boundary_centroid_x, right%boundary_centroid_x)) .and. &
      all(checkpoint_real_matches(left%boundary_centroid_y, right%boundary_centroid_y)) .and. &
      all(checkpoint_real_matches(left%boundary_centroid_z, right%boundary_centroid_z)) .and. &
      all(checkpoint_real_matches(left%boundary_normal_x, right%boundary_normal_x)) .and. &
      all(checkpoint_real_matches(left%boundary_normal_y, right%boundary_normal_y)) .and. &
      all(checkpoint_real_matches(left%boundary_normal_z, right%boundary_normal_z)) .and. &
      all(checkpoint_real_matches(left%boundary_normal_integral_x, &
        right%boundary_normal_integral_x)) .and. &
      all(checkpoint_real_matches(left%boundary_normal_integral_y, &
        right%boundary_normal_integral_y)) .and. &
      all(checkpoint_real_matches(left%boundary_normal_integral_z, &
        right%boundary_normal_integral_z)) .and. &
      all(left%cell_type == right%cell_type)
  end function geometries_match

  subroutine allocate_geometry_storage(geometry, status)
    type(eb_geometry_3d), intent(inout) :: geometry
    integer, intent(out) :: status

    status = 1
    if (geometry%nx < 1 .or. geometry%ny < 1 .or. geometry%nz < 1) return
    allocate(geometry%volume_fraction(geometry%nx, geometry%ny, geometry%nz), &
      geometry%cell_centroid_x(geometry%nx, geometry%ny, geometry%nz), &
      geometry%cell_centroid_y(geometry%nx, geometry%ny, geometry%nz), &
      geometry%cell_centroid_z(geometry%nx, geometry%ny, geometry%nz), &
      geometry%x_face_fraction(0:geometry%nx, geometry%ny, geometry%nz), &
      geometry%y_face_fraction(geometry%nx, 0:geometry%ny, geometry%nz), &
      geometry%z_face_fraction(geometry%nx, geometry%ny, 0:geometry%nz), &
      geometry%x_face_centroid_y(0:geometry%nx, geometry%ny, geometry%nz), &
      geometry%x_face_centroid_z(0:geometry%nx, geometry%ny, geometry%nz), &
      geometry%y_face_centroid_x(geometry%nx, 0:geometry%ny, geometry%nz), &
      geometry%y_face_centroid_z(geometry%nx, 0:geometry%ny, geometry%nz), &
      geometry%z_face_centroid_x(geometry%nx, geometry%ny, 0:geometry%nz), &
      geometry%z_face_centroid_y(geometry%nx, geometry%ny, 0:geometry%nz), &
      geometry%boundary_area(geometry%nx, geometry%ny, geometry%nz), &
      geometry%boundary_centroid_x(geometry%nx, geometry%ny, geometry%nz), &
      geometry%boundary_centroid_y(geometry%nx, geometry%ny, geometry%nz), &
      geometry%boundary_centroid_z(geometry%nx, geometry%ny, geometry%nz), &
      geometry%boundary_normal_x(geometry%nx, geometry%ny, geometry%nz), &
      geometry%boundary_normal_y(geometry%nx, geometry%ny, geometry%nz), &
      geometry%boundary_normal_z(geometry%nx, geometry%ny, geometry%nz), &
      geometry%boundary_normal_integral_x(geometry%nx, geometry%ny, geometry%nz), &
      geometry%boundary_normal_integral_y(geometry%nx, geometry%ny, geometry%nz), &
      geometry%boundary_normal_integral_z(geometry%nx, geometry%ny, geometry%nz), &
      geometry%cell_type(geometry%nx, geometry%ny, geometry%nz), stat=status)
  end subroutine allocate_geometry_storage

  subroutine write_real_array_3d(unit, values, status)
    integer, intent(in) :: unit
    real(dp), intent(in) :: values(:, :, :)
    integer, intent(inout) :: status
    integer :: i, j, k

    if (status /= 0) return
    do k = lbound(values, 3), ubound(values, 3)
      do j = lbound(values, 2), ubound(values, 2)
        do i = lbound(values, 1), ubound(values, 1)
          write(unit, '(es27.18e3)', iostat=status) values(i, j, k)
          if (status /= 0) return
        end do
      end do
    end do
  end subroutine write_real_array_3d

  subroutine read_real_array_3d(unit, values, status)
    integer, intent(in) :: unit
    real(dp), intent(out) :: values(:, :, :)
    integer, intent(inout) :: status
    integer :: i, j, k

    if (status /= 0) return
    do k = lbound(values, 3), ubound(values, 3)
      do j = lbound(values, 2), ubound(values, 2)
        do i = lbound(values, 1), ubound(values, 1)
          read(unit, *, iostat=status) values(i, j, k)
          if (status /= 0) return
        end do
      end do
    end do
  end subroutine read_real_array_3d

  subroutine write_integer_array_3d(unit, values, status)
    integer, intent(in) :: unit
    integer, intent(in) :: values(:, :, :)
    integer, intent(inout) :: status
    integer :: i, j, k

    if (status /= 0) return
    do k = lbound(values, 3), ubound(values, 3)
      do j = lbound(values, 2), ubound(values, 2)
        do i = lbound(values, 1), ubound(values, 1)
          write(unit, '(i0)', iostat=status) values(i, j, k)
          if (status /= 0) return
        end do
      end do
    end do
  end subroutine write_integer_array_3d

  subroutine read_integer_array_3d(unit, values, status)
    integer, intent(in) :: unit
    integer, intent(out) :: values(:, :, :)
    integer, intent(inout) :: status
    integer :: i, j, k

    if (status /= 0) return
    do k = lbound(values, 3), ubound(values, 3)
      do j = lbound(values, 2), ubound(values, 2)
        do i = lbound(values, 1), ubound(values, 1)
          read(unit, *, iostat=status) values(i, j, k)
          if (status /= 0) return
        end do
      end do
    end do
  end subroutine read_integer_array_3d

  pure elemental logical function checkpoint_real_matches(stored, expected) &
      result(matches)
    real(dp), intent(in) :: stored, expected
    real(dp) :: tolerance

    tolerance = 64.0_dp * epsilon(1.0_dp) * &
      max(1.0_dp, abs(stored), abs(expected))
    matches = ieee_is_finite(stored) .and. ieee_is_finite(expected) .and. &
      abs(stored - expected) <= tolerance
  end function checkpoint_real_matches

end module reactive_eb_3d_checkpoint_mod
